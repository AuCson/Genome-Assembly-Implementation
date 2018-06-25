"""
evalutate N50
"""
import argparse
from gene import read_fasta

def N50(filename):
    l = read_fasta(filename)
    l.sort(key=lambda x: -len(x))
    sum_len = sum([len(_) for _ in l])
    s = 0
    for contig in l:
        s += len(contig)
        if s > sum_len / 2:
            return len(contig)
    return -1

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    args = parser.parse_args()
    score = N50(args.file)
    print(score)

if __name__ == '__main__':
    main()
        