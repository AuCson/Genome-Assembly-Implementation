import argparse

def filter_line(l, thres=100):
    i = 0
    ret = []
    while i < len(l):
        header = l[i]
        seq = l[i+1]
        if len(seq.strip()) >= thres:
            ret.append(header.strip()+'\n')
            ret.append(seq.strip()+'\n')
        i += 2
    return ret


def cutline(fa, fo):
    l = fa.readlines()
    ret = filter_line(l,  100)
    fo.writelines(ret)

def recur_merge(file_list):
    ret = []
    max_len = 0
    lim = [6000, 7600, 7600, 7600]
    for i,file in enumerate(file_list):
        lines = open(file).readlines()
        thres = min(max_len-1000,9200)
        l = filter_line(lines, thres)
        ret.extend(l)
        if l:
            ml = max([len(_) for _ in l])
            if max_len < ml:
                max_len = ml
            print(max_len)
    return ret


if __name__ == '__main__':
    merge = ['full1_contig17.txt','full1_good.txt','full1.txt','full1_contig51.txt','full1_contig53.txt','full1_contig55.txt',
    'full1_contig57.txt','full1_contig59.txt','full1_contig61.txt',
    'full1_contig63.txt','full1_contig65.txt','full1_contig67.txt','full1_contig69.txt','full1_contig71.txt',
    'full1_contig73.txt','full1_contig75.txt','full1_contig77.txt']
    merge = ['result/' + _ for _ in merge]
    line = recur_merge(merge)
    f = open('result/full1_ultramerge.txt','w')
    f.writelines(line)