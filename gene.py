CHAR_MAP = {
    'A':'T',
    'T':'A',
    'C':'G',
    'G':'C',
    'X':'X'
}

def compliment(s):
    v = []
    for c in s:
        v.append(CHAR_MAP.get(c,''))
    return ''.join(c)

def rev_complement(s):
    return ''.join(reversed(compliment(s)))

def canon(s):
    cmpl = compliment(s)
    return s if s < cmpl else cmpl

def compif(s, i):
    if i != '1' and i != '2':
        raise ValueError()
    return s if i == '1' else rev_complement(s)