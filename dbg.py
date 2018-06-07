"""

Author: Xisen Jin
Date: 2018.05.27

De Brujin Graph

"""

from gene import *
import collections

class Edge:
    def __init__(self, seq1, seq2, ori):
        # seq1 --> seq2
        self.seq1 = seq1
        self.seq2 = seq2
        self.ori = ori # '11','12','21','22'
        self.cov = 1

class Vertex:
    def __init__(self, seq):
        self.seq = seq
        self.cov = 1
        self.out_edge = {}

    def rev(self):
        return rev_complement(self.seq)

class DBG:
    """
    Data structure
    
    seq --O(1)--> vertex
    vertex + seq --O(1) --> edge
    
    """
    def __init__(self, k):
        self.k = k
        self.v = {}

    def add_vertex(self, kmer):
        kmer = canon(kmer) # store the smaller one
        if kmer in self.v:
            self.v[kmer].cov += 1
        else:
            self.v[kmer] = Vertex(kmer)

    def add_edge(self, kmer1, kmer2):
        """
        edge: kmer1 to kmer2
        """
        ckmer1 = canon(kmer1)
        ckmer2 = canon(kmer2)
        ori1 = '1' if ckmer1 == kmer1 else '2'
        ori2 = '1' if ckmer2 == kmer2 else '2'
        if ckmer2 not in self.v[ckmer1].out_edge:
            self.v[ckmer1].out_edge[ckmer2] = Edge(ckmer1, ckmer2, ori1 + ori2)
        else:
            self.v[ckmer1].out_edge[ckmer2].cov += 1
        if ckmer1 not in self.v[ckmer2].out_edge:
            self.v[ckmer2].out_edge[ckmer1] = Edge(ckmer2, ckmer1, revori(ori1 + ori2))
        else:
            self.v[ckmer2].out_edge[ckmer1].cov += 1

    def add_read(self, read):
        if len(read) <= self.k:
            return
        kmer1, kmer2 = None, None
        for i in range(len(read)-self.k):
            kmer1 = read[i:i+self.k]
            kmer2 = read[i+1:i+self.k+1]
            self.add_vertex(kmer1)
            self.add_vertex(kmer2)
            self.add_edge(kmer1, kmer2)

    def dfs_wrapper(self, start_node, visit, start_ori):
        #path = list(start_node.seq[:-1])
        full_path = list(start_node.seq[:-1])
        node_stack = collections.deque()
        node_stack.append((start_node, start_ori))
        while len(node_stack):
            node, ori = node_stack.pop()
            full_path.extend(list(node.seq[self.k-1:]))
            visit[node.seq] = True
            for edge in node.out_edge.values():
                if edge.ori[0] == ori and edge.seq2 not in visit:
                    node_stack.append((self.v[edge.seq2], edge.ori[1]))
        return full_path

    def get_contig_wrapper(self, start_node, visit, start_ori):

        def is_self_ring(v):
            return v.seq in v.out_edge

        node_stack = collections.deque()
        node_stack.append((start_node, start_ori))
        contigs = []
        a = 0
        # stack stores nodes whose degree != 2
        while len(node_stack):
            a += 1
            b = 0
            node, ori = node_stack.pop()
            visit.add((node.seq, ori))
            #print(len(self.v), len(visit), a)
            for edge in node.out_edge.values():
                # travel until degree > 2 or degree = 1

                if edge.ori[0] == ori:
                    path = []
                    t_ori = edge.ori[1]
                    v = self.v[edge.seq2]
                    while True:
                        b += 1
                        #print(len(self.v), len(visit), a,b)
                        path.append(compif(v.seq, t_ori)[-1])
                        if (v.seq, t_ori) in visit:
                            break
                        visit.add((v.seq, t_ori))
                        if len(v.out_edge) != 2 or is_self_ring(v):
                            if len(v.out_edge) >= 2:
                                node_stack.append((v, t_ori))
                            break
                        else:
                            next_v = None
                            for edge_ in v.out_edge.values():
                                if edge_.ori[0] == t_ori:
                                    next_v = self.v[edge_.seq2]
                                    t_ori = edge_.ori[1]
                                    break
                            if next_v is None:
                                break
                            v = next_v
                    contigs.append(compif(node.seq, ori) + ''.join(path))
        return contigs           

    def dfs_graph(self):
        visit = {}
        all_path = []
        while True:
            start_node = None
            for k in self.v:
                if k not in visit:
                    start_node = self.v[k]
                    break
            if not start_node:
                break
            path = self.dfs_wrapper(start_node, visit, '1')
            all_path.append(''.join(path)+'\n')
        return all_path
        

    def get_contig_graph(self):
        """
        dfs and store any unambiguous contig
        """
        visit = set()
        all_path = []
        c = 0
        while True:
            c += 1
            print(c,len(self.v),len(visit))
            if (len(visit) >= len(self.v)):
                break
            start_node = None
            for k in self.v:
                if (k,'1') not in visit and len(self.v[k].out_edge)!=2:
                    start_node = self.v[k]
                    break
            if not start_node:
                break
            paths = self.get_contig_wrapper(start_node, visit, '1')
            all_path.extend(paths)

            for k in self.v:
                if (k,'2') not in visit and len(self.v[k].out_edge)!=2:
                    start_node = self.v[k]
                    break
            if not start_node:
                break
            paths = self.get_contig_wrapper(start_node, visit, '2')
            all_path.extend(paths)
        return all_path

def write_fa(f, lines):
    for i,line in enumerate(lines):
        #if len(line) < 100:
        #    continue
        f.write('>{} length {} xxx\n'.format(i*2+1, len(line.strip())))
        #f.write('>xxx\n')
        f.write(line.strip()+'\n')
        #f.write('>xxx\n')
        #f.write(rev_complement(line.strip())+'\n')

if __name__ == '__main__':
    g = DBG(k=25)
    reads = read_fasta('data/data2/short_1.fasta')
    reads += read_fasta('data/data2/short_2.fasta')
    #reads = ['ATCGA']
    for i,read in enumerate(reads):
        g.add_read(read)
        print(i)
        #if i == 10000:
        #    break
    contigs = g.get_contig_graph()
    #contigs = g.dfs_graph()
    f = open('result/data2.txt','w')
    write_fa(f, contigs)

'''
    def graph_simplification(self):
        """
        simplify the graph by merging nodes
        implemented with BFS
        """
        contigs = []
        visit = {} # key is seq, not vertex
        print(len(self.v))
        def simplify_path(v, in_ori):
            """
            simplify a chain that starts from the v
            return a new vertex and final visited node, 
            """
            turn_visit = {}
            seqs = []
            final_edge = None
            while len(v.out_edge) <= 2:
                if not seqs:
                    seqs.extend(list(compif(v.seq, in_ori)))
                else:
                    seqs.append(compif(v.seq, in_ori)[-1])
                visit[v.seq] = True
                turn_visit[v.seq] =
                next_edge = None
                for edge in v.out_edge.values():
                    if in_ori == edge.ori[0] and edge.seq2 not in turn_visit: # compatible
                        next_edge = edge
                        break
                if next_edge is None:
                    v = None
                    final_edge = None
                    break
                final_edge = v.seq
                v = self.v[next_edge.seq2]
                in_ori = edge.ori[1]
                assert(final_edge in v.out_edge)
            
            return v, ''.join(seqs), final_edge

        def merge(s,e,new_seq, first_edge_seq,final_edge_seq):
            self.add_vertex(new_seq)
            visit[canon(new_seq)] = True
            s.out_edge.pop(first_edge_seq)
            self.add_edge(s.seq, new_seq)
            if e is not None:
                e.out_edge.pop(final_edge_seq)
                self.add_edge(new_seq, e.seq)

        cnt = 0
        while True:
            # initialize
            cnt += 1
            q = collections.deque()
            first_item = None
            for k in self.v:
                if len(self.v[k].out_edge) != 2 and k not in visit:
                    first_item = self.v[k]
                    break
            if first_item is None:
                break
            q.append(first_item)
            cnt2 = 0
            while len(q):
                print('loop %d %d %d %d' % (cnt, cnt2,len(self.v), len(visit)))
                cnt2 += 1
                v = q.popleft()
                if v.seq in visit:
                    continue
                visit[v.seq] = True
                for edge in v.out_edge.values():
                    ori = edge.ori[1]
                    next_v, seq, final_edge_seq = simplify_path(self.v[edge.seq2], ori)
                    seq_head = compif(edge.seq1, edge.ori[0])[0]
                    seq_tail = ''
                    if seq:
                        merge(v, next_v, seq, edge.seq2, final_edge_seq)
                    if next_v is not None and next_v.seq not in visit:
                        q.append(next_v)
                    if next_v is not None:
                        for edge in next_v.out_edge.values():
                            if edge.seq2 == final_edge_seq:
                                seq_tail = compif(edge.seq1, '1' if edge.ori[0] == '2' else '1')[-1]
                    contig = seq_head + seq + seq_tail
                    contigs.append(contig)
        return contigs
'''
        
    

        









