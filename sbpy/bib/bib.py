from collections import OrderedDict

__all__ = ['Bib']

class Bib():
    """Bibliography class

    References for specific tasks are provided as bibcode elements that can be queried at `ADS`_

    .. _ADS: http://adsabs.harvard.edu
    """

    def __init__(self):
        self.bib = OrderedDict()

    def __setitem__(self, task, payload):
        """Set references for a specific task

        Parameters
        ----------
        task : str, mandatory
            name of the task generating the reference
        payload : dict, mandatory
            a dictionary with bibcodes for several aspects of the task, e.g., 
            one could cite one paper for the general `method` and one for the
            actual `implementation` used in SBPy
        """
        self.bib[task] = payload

    def __getitem__(self, ident):
        """Return the references for a specific task or for the n-th task""" 
        if isinstance(ident, str):
            return self.bib[ident]
        elif isinstance(ident, int):
            key = list(self.bib.keys())[idx]
            return {'task': key, 'references': self.bib[key]}
        
    def __iter__(self):
        return self.bib
            
    def __len__(self):
        return len(self.bib)
    
    @property
    def tasks(self):
        return list(self.bib.keys())
    
    def to_text(self):
        """convert bibcodes to human readable text

        not yet implemented
        """
        output = ''
        for ref in self.bib: #__iter__():
            output += '{:s}:\n'.format(ref)
            for key, val in self.__getitem__(ref).items():
                # transform val to ads(val)
                output += '  {:s}: {:s}\n'.format(key, val)
        return output

    def to_bibtex(self):
        """ convert bibcodes to LATeX bibtex

        not yet implemented
        """
        output = ''
        for ref in self.bib: #__iter__():
            for key, val in self.__getitem__(ref).items():
                # transform val to ads.bibtex(val)
                output += '% {:s}/{:s}:\n{:s}\n'.format(ref, key, val)
        return output


