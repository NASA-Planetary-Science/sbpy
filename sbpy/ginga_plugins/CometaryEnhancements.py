# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Skeleton example of a Ginga local plugin called 'CometaryEnhancements'.

With sbpy installed, this plugin should be automatically discovered by
Ginga and available in the Operations menu.

"""

from warnings import warn
from sbpy.exceptions import OptionalPackageUnavailable

try:
    from ginga.GingaPlugin import LocalPlugin
    from ginga.gw import Widgets
except ImportError:
    warn('ginga is not present: CometaryEnhancements will not run.',
         OptionalPackageUnavailable)

    class LocalPlugin:
        pass
    Widgets = None

try:
    import photutils
except ImportError:
    warn('photutils is not present: CometaryEnhancements centroiding disabled.',
         OptionalPackageUnavailable)
    photutils = None

from sbpy.imageanalysis import CometaryEnhancement


class CometaryEnhancements(LocalPlugin):
    """Ginga plugin for on-the-fly cometary image enhancements."""

    def __init__(self, fv, fitsimage):
        """
        This method is called when the plugin is loaded for the  first
        time.  ``fv`` is a reference to the Ginga (reference viewer) shell
        and ``fitsimage`` is a reference to the specific ImageViewCanvas
        object associated with the channel on which the plugin is being
        invoked.
        You need to call the superclass initializer and then do any local
        initialization.
        """
        super(CometaryEnhancements, self).__init__(fv, fitsimage)

        # your local state and initialization code goes here
        self.enhancement_options = ['1/rho']

    def build_gui(self, container):
        """
        This method is called when the plugin is invoked.  It builds the
        GUI used by the plugin into the widget layout passed as
        ``container``.
        This method may be called many times as the plugin is opened and
        closed for modal operations.  The method may be omitted if there
        is no GUI for the plugin.

        This specific example uses the GUI widget set agnostic wrappers
        to build the GUI, but you can also just as easily use explicit
        toolkit calls here if you only want to support one widget set.
        """
        top = Widgets.VBox()
        top.set_border_width(4)

        # this is a little trick for making plugins that work either in
        # a vertical or horizontal orientation.  It returns a box container,
        # a scroll widget and an orientation ('vertical', 'horizontal')
        vbox, sw, orientation = Widgets.get_oriented_box(container)
        vbox.set_border_width(4)
        vbox.set_spacing(2)

        # Take a text widget to show some instructions
        self.msg_font = self.fv.get_font("sansFont", 12)
        tw = Widgets.TextArea(wrap=True, editable=False)
        tw.set_font(self.msg_font)
        self.tw = tw

        # Frame for instructions and add the text widget with another
        # blank widget to stretch as needed to fill emp
        frame = Widgets.Expander("Instructions")
        frame.set_widget(tw)
        vbox.add_widget(frame, stretch=0)

        frame = Widgets.Frame('Center')
        hbox = Widgets.HBox()
        w, b = Widgets.build_info(
            (('X center:', 'label', 'X center', 'entry'),
             ('Y center:', 'label', 'Y center', 'entry'),
             ('Background:', 'label', 'background', 'entry'),)
        )
        b.background.set_text('0.0')
        self.w.update(b)
        hbox.add_widget(w)

        w, b = Widgets.build_info(
            (('Use FOV center', 'button'),
             ('Centroid', 'button'),
             ('Centroid box:', 'label', 'centroid box', 'entry'),)
        )
        b.use_fov_center.add_callback('activated', self.use_fov_center_cb)
        # centroiding requires photutils
        if photutils:
            b.centroid.add_callback('activated', self.centroid_cb)
            b.centroid_box.set_text('11')
        else:
            b.centroid.set_enabled(False)
            b.centroid_box.set_enabled(False)

        self.w.update(b)
        hbox.add_widget(w)

        frame.set_widget(hbox)
        vbox.add_widget(frame)

        # Enhancement tabs
        tabw = Widgets.TabWidget()
        self.w.tabw = tabw

        # 1/rho
        rho_vbox = Widgets.VBox()
        widgets = (('No options', 'label'),
                   ('Enhance', 'button'))
        w, b = Widgets.build_info(widgets)
        b.enhance.add_callback('activated', self.rho_cb)
        self.w.rho = b
        rho_vbox.add_widget(w)
        tabw.add_widget(rho_vbox, title='1/rho')

        vbox.add_widget(tabw)

        # scroll bars will allow lots of content to be accessed
        top.add_widget(sw, stretch=1)

        # A button box that is always visible at the bottom
        btns = Widgets.HBox()
        btns.set_spacing(3)

        # Add a close button for the convenience of the user
        btn = Widgets.Button("Close")
        btn.add_callback('activated', lambda w: self.close())
        btns.add_widget(btn, stretch=0)
        btns.add_widget(Widgets.Label(''), stretch=1)
        top.add_widget(btns, stretch=0)

        # Add our GUI to the container
        container.add_widget(top, stretch=1)
        # NOTE: if you are building a GUI using a specific widget toolkit
        # (e.g. Qt) GUI calls, you need to extract the widget or layout
        # from the non-toolkit specific container wrapper and call on that
        # to pack your widget, e.g.:
        # cw = container.get_widget()
        # cw.addWidget(widget, stretch=1)

    def rho_cb(self, w):
        try:
            xc = float(self.w.x_center.get_text())
            yc = float(self.w.y_center.get_text())
        except ValueError:
            return

        original = self.fitsimage.get_image()
        new_image = original.copy()

        try:
            bg = float(self.w.background.get_text())
        except ValueError:
            bg = 0

        enhancer = CometaryEnhancement(new_image.get_data() - bg, (yc, xc))
        new_image.set_data(enhancer.rho_norm())

        chname = self.fv.get_current_channel().name + '(1/rho)'
        self.fv.add_image('1/rho enhanced', new_image, chname=chname)

    def centroid_cb(self, w):
        im = self.fitsimage.get_image().get_data()

        try:
            xc = int(float(self.w.x_center.get_text()))
            yc = int(float(self.w.y_center.get_text()))
            box = int(float(self.w.centroid_box.get_text()))
        except ValueError:
            return

        if box < 3:
            box = 3
            self.w.centroid_box.set_text('3')
        elif box > min(xc, yc, im.shape[0] - yc, im.shape[1] - xc):
            box = int(min(xc, yc, im.shape[0] - yc, im.shape[1] - xc))
            self.w.centroid_box.set_text(str(box))

        x0 = xc - box // 2
        y0 = yc - box // 2

        subim = im[y0:y0 + box + 1, x0:x0 + box + 1]
        cxy = photutils.centroid_2dg(subim)
        self.w.x_center.set_text('{:.2f}'.format(cxy[0] + x0))
        self.w.y_center.set_text('{:.2f}'.format(cxy[1] + y0))

    def use_fov_center_cb(self, w):
        """Use the field of view center for the enhancement."""
        xy = self.fv.get_viewer(self.chname).get_pan()
        self.w.x_center.set_text('{:.2f}'.format(xy[0]))
        self.w.y_center.set_text('{:.2f}'.format(xy[1]))

    def close(self):
        """
        Example close method.  You can use this method and attach it as a
        callback to a button that you place in your GUI to close the plugin
        as a convenience to the user.
        """
        self.fv.stop_local_plugin(self.chname, str(self))
        return True

    def start(self):
        """
        This method is called just after ``build_gui()`` when the plugin
        is invoked.  This method may be called many times as the plugin is
        opened and closed for modal operations.  This method may be omitted
        in many cases.
        """
        self.tw.set_text("""Select target center and choose enhancement.""")
        self.resume()

    def pause(self):
        """
        This method is called when the plugin loses focus.
        It should take any actions necessary to stop handling user
        interaction events that were initiated in ``start()`` or
        ``resume()``.
        This method may be called many times as the plugin is focused
        or defocused.  It may be omitted if there is no user event handling
        to disable.
        """
        pass

    def resume(self):
        """
        This method is called when the plugin gets focus.
        It should take any actions necessary to start handling user
        interaction events for the operations that it does.
        This method may be called many times as the plugin is focused or
        defocused.  The method may be omitted if there is no user event
        handling to enable.
        """
        pass

    def stop(self):
        """
        This method is called when the plugin is stopped.
        It should perform any special clean up necessary to terminate
        the operation.  The GUI will be destroyed by the plugin manager
        so there is no need for the stop method to do that.
        This method may be called many times as the plugin is opened and
        closed for modal operations, and may be omitted if there is no
        special cleanup required when stopping.
        """
        pass

    def redo(self):
        """
        This method is called when the plugin is active and a new
        image is loaded into the associated channel.  It can optionally
        redo the current operation on the new image.  This method may be
        called many times as new images are loaded while the plugin is
        active.  This method may be omitted.
        """
        pass

    def __str__(self):
        """
        This method should be provided and should return the lower case
        name of the plugin.
        """
        return 'cometaryenhancements'
