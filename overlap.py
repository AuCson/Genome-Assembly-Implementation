import numpy as np
import time

def find_back_dp(path,lp,ls):
    i,j = lp,ls
    p = []
    while i >= 0 and j >= 0:
        i,j,sig = backward(path, i, j)
        p.append((i,j))
        if sig == 0:
            break
    p = [_ for _ in reversed(p)]
    print(p)
    return i,j

def backward(path, i, j):
    d = {
        0: (i,j,0),
        1: (i-1, j-1, 1),
        2: (i-1, j, 1),
        3: (i, j-1, 1)
    }
    return d[path[i,j]]

def approximate_find(pat, seq, TOL=2):
    """
    find an approximate overlapping sequence. 
    """
    lp, ls = len(pat), len(seq)
    dp = np.zeros((lp+1, ls+1), dtype=np.int32)
    max_dp = -1
    max_idx = 0,0
    toler = np.full((lp+1, ls+1), TOL)
    path = np.zeros((lp+1, ls+1), dtype=np.int32) # 1 left-up; 2 left; 3 up;
    i,j = 1,1
    while i <= lp:
        j = 1
        while j<= ls:

            if pat[i-1] == seq[j-1]:
                dp[i,j] = dp[i-1,j-1] + 1
                path[i,j] = 1
                toler[i,j] = toler[i-1,j-1] + 0.25
                if dp[i,j] > max_dp:
                    max_idx = i,j
                    max_dp = dp[i,j]
            if toler[i-1,j] > 0:
                if dp[i-1,j] > dp[i,j]:
                    dp[i,j] = dp[i-1,j]
                    path[i,j] = 2
                    toler[i,j] = toler[i-1,j] - 1

            if toler[i,j-1] > 0:
                if dp[i,j-1] > dp[i,j]:
                    dp[i,j] = dp[i,j-1]
                    path[i,j] = 3
                    toler[i,j] = toler[i,j-1] - 1

            j += 1
        i += 1

    p = find_back_dp(path, max_idx[0], max_idx[1])
    print(dp.max())

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

def main():
    fa = read_fasta('./data/data1/graph_prefix.contig')
    fb = read_fasta('./data/data1/long.fasta')
    print(len(fa[-1]),len(fb[0]))
    a = time.time()
    approximate_find(fa[-1],fb[0])
    print(time.time() - a)

if __name__ == '__main__':
    main()
    

            