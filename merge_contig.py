import argparse

def filter_line(l, thres=100):
    i = 0
    ret = []
    while i < len(l):
        header = l[i]
        seq = l[i+1]
        if len(seq.strip()) <= thres:
            ret.append(header)
            ret.append(seq)
        i += 2
    return ret

def concat(l, fo):
    ret = []
    seq = ''
    i = 0
    while i < len(l):
        seq += l[i+1].strip()
        if len(seq) > 6000:
            ret.append('> %d\n' % len(seq))
            ret.append(seq + '\n') 
            seq = ''
        i += 2
    if seq:
        ret.append('> %d\n' % len(seq))
        ret.append(seq +'\n' )
    fo.writelines(ret)
        

def merge(fa, fb, fo):
    l = fa.readlines()
    la = filter_line(l, thres=1000000)
    lb = filter_line(fb.readlines(), thres=10)
    #assert(la == l and lb == [] and 1)

    ret = la + lb
    fo.writelines(ret)

if __name__ == '__main__':
    f1 = open('result/data_soap2.txt')
    f2 = open('result/data2_contig.txt')
    fo = open('result/data2_merge.txt','w')
    merge(f1, f2, fo)
    #concat(f1.readlines(), fo)
