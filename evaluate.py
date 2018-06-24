"""
evalutate N50
"""

def N50(filename):
    f = open(filename)
    l = f.readlines()
    l.sort(key=lambda x: -len(x))
    if l == 0:
        return True