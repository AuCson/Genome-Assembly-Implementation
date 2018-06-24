import argparse

def filter_line(l, thres=100, chunk=0):
    i = 0
    ret = []
    while i < len(l):
        header = l[i]
        seq = l[i+1]
        if len(seq.strip()) >= thres:
            if not chunk:
                ret.append(header)
                ret.append(seq)
            else:
                for j in range(0, len(seq),chunk- 1):
                    ret.append(header)
                    ret.append(seq[j:min(j+chunk, len(seq))].strip()+'\n')
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
        

def cat(fa, fb, fo):
    l = fa.readlines()
    la = filter_line(l, thres=0)
    lb = filter_line(fb.readlines(), thres=0, chunk=50)
    #assert(la == l and lb == [] and 1)

    ret = la + lb
    fo.writelines(ret)

def cutline(fa, fo):
    l = fa.readlines()
    ret = filter_line(l, 0, 100)
    fo.writelines(ret)

def extend(fa, fo):
    for l in fa.readlines():
        if len(l) > 1000:
            l = l.strip() + ''.join(['A'] * 500) + '\n'
            fo.write(l)
        else:
            fo.write(l)

if __name__ == '__main__':
    #f1 = open('result/data_soap3.txt')
    f2 = open('result/full4_contig.txt')
    fo = open('result/data4_cat.txt','w')
    #merge(f1, f2, fo)
    #cutline(f2, fo)
    extend(f2, fo)
