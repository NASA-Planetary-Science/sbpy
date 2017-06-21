# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""Skeleton example of a Ginga local plugin called 'CometaryEnhancements'.

With sbpy installed, this plugin should be automatically discovered by
Ginga and available in the Operations menu.

"""

from ginga import GingaPlugin
from ginga.gw import Widgets

class CometaryEnhancements(GingaPlugin.LocalPlugin):

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
        self.image = None

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
             ('Y center:', 'label', 'Y center', 'entry'))
        )
        self.w.update(b)
        hbox.add_widget(w)

        w, b = Widgets.build_info(
            (('Pick center', 'button'),
             ('Centroid', 'button'))
        )
        self.w.update(b)
        hbox.add_widget(w)

        frame.set_widget(hbox)
        vbox.add_widget(frame)
        
        captions = (('Enhancement:', 'label', 'Enhancement', 'combobox',
                     'Enhance', 'button'),)
        w, b = Widgets.build_info(captions)
        self.w.update(b)
        b.enhance.add_callback('activated', self.enhance_cb)

        for name in self.enhancement_options:
            b.enhancement.append_text(name)
        b.enhancement.set_tooltip('Cometary enhancement method')
        b.enhancement.set_index(0)
        b.enhancement.add_callback('activated', lambda w, i: self.set_enhancement_cb())

        vbox.add_widget(w)

        frame = Widgets.Frame('Enhancement options')
        self.w.enhancementsvbox = Widgets.VBox()
        frame.set_widget(self.w.enhancementsvbox, stretch=1)
        vbox.add_widget(frame)
        self.set_enhancement_cb()
        
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
        #cw = container.get_widget()
        #cw.addWidget(widget, stretch=1)

    def enhance_cb(self, w):
        import numpy as np
        
        try:
            xc = float(self.w.x_center.get_text())
            yc = float(self.w.y_center.get_text())
        except ValueError:
            return

        if self.image is None:
            self.image = self.fitsimage.get_image()

        y, x = np.indices(self.image.shape, float)
        rho = np.sqrt((x - xc)**2 + (y - yc)**2)
        
        enhanced = self.image.copy()
        enhanced.set_data(enhanced.get_data() * rho)
        chname = self.fv.get_current_channel().name + '(1/rho)'
        self.fv.add_image('1/rho enhanced', enhanced, chname=chname)
        
    def set_enhancement_cb(self):
        i = self.w.enhancement.get_index()

        # remove old parameters
        self.w.enhancementsvbox.remove_all()

        # create new parameters
        captions = (('None', 'label'),)
        w, b = Widgets.build_info(captions)
        self.w.enhancementsvbox.add_widget(w)

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
        This method may be called many  times as the plugin is opened and
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
