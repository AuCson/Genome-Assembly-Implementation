CHAR_MAP = {
    'A':'T',
    'T':'A',
    'C':'G',
    'G':'C',
    'X':'X'
}

def complement(s):
    v = [None] * len(s)
    for i,c in enumerate(s):
        if c == 'A': v[i] = 'T'
        elif c == 'T': v[i] = 'A'
        elif c == 'G': v[i] = 'C'
        elif c == 'C': v[i] = 'G'
    return ''.join(v)

def rev_complement(s):
    return ''.join(reversed(complement(s)))

def canon(s):
    cmpl = rev_complement(s)
    return s if s < cmpl else cmpl

def compif(s, i):
    if i != '1' and i != '2':
        raise ValueError()
    return s if i == '1' else rev_complement(s)

def revori(ori):
    ORI_MAP = {
        '11':'22',
        '22':'11',
        '12':'12',
        '21':'21',
        '1':'2',
        '2':'1'
    }
    return ORI_MAP[ori]

def read_fasta(path):
    f = open(path)
    l = []
    s = ''
    for line in f.readlines():
        if line.startswith('>'):
            if s:
                l.append(s)
            s = ''
        else:
            s += line.strip()
    if s:
        l.append(s)
    return l
