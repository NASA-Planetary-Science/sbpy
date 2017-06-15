from collections import OrderedDict

class BibList():
    def __init__(self):
        self.bib = OrderedDict()

    def __setitem__(self, task, payload):
        self.bib[task] = payload

    def __getitem__(self, ident):
        if isinstance(ident, str):
            return self.bib[ident]
        elif isinstance(ident, int):
            key = list(self.bib.keys())[idx]
            return {'task': key, 'references': self.bib[key]}
        
    def __iter__(self):
        return self.bib
        # for ref in list(self.bib):
        #     yield ref
            
    def __len__(self):
        return len(self.bib)
    
    # def __repr__(self):
    #     return str(self.bib)

    @property
    def tasks(self):
        return list(self.bib.keys())
    
    def to_text(self):
        output = ''
        for ref in self.bib: #__iter__():
            output += '{:s}:\n'.format(ref)
            for key, val in self.__getitem__(ref).items():
                # transform val to ads(val)
                output += '  {:s}: {:s}\n'.format(key, val)
        return output

    def to_bibtex(self):
        output = ''
        for ref in self.bib: #__iter__():
            for key, val in self.__getitem__(ref).items():
                # transform val to ads.bibtex(val)
                output += '% {:s}/{:s}:\n{:s}\n'.format(ref, key, val)
        return output


# example/test code
# -----------------

# from sbpy import Reference

# test = Reference.BibList()
# test['stuff'] = {'method': 'A', 'implementation': 'B'}
# test['extrastuff'] = {'method': 'C', 'implementation': 'D'}

# print(len(test))

# print(test.tasks)

# print(test.to_text())

# print(test.to_bibtex())
