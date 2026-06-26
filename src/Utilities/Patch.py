from DataBase import *

def nodelist_getattr(self, name):
    try:
        return self.getFieldDouble(name)
    except:
        pass
    try:
        return self.getFieldComplex(name)
    except:
        pass
    for dim in [1, 2, 3]:
        try:
            return getattr(self, f"getFieldVector{dim}d")(name)
        except:
            pass
    raise AttributeError(f"Field '{name}' not found.")

# Patch the method onto the generated class
NodeList.__getattr__ = nodelist_getattr