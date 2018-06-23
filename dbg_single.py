"""

Author: Xisen Jin
Date: 2018.05.27

De Brujin Graph

"""

from gene import *
import collections

class Edge:
    def __init__(self, seq1, seq2):
        # seq1 --> seq2
        self.seq1 = seq1
        self.seq2 = seq2
        self.cov = 1

class Vertex:
    def __init__(self, seq):
        self.seq = seq
        self.cov = 1
        self.out_edge = {}

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
        if kmer in self.v:
            self.v[kmer].cov += 1
        else:
            self.v[kmer] = Vertex(kmer)

    def add_edge(self, kmer1, kmer2):
        """
        edge: kmer1 to kmer2
        """
        if kmer2 not in self.v[kmer1].out_edge:
            self.v[kmer1].out_edge[kmer2] = Edge(kmer1, kmer2)
        else:
            self.v[kmer1].out_edge[kmer2].cov += 1

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

            self.add_vertex(rev_complement(kmer1))
            self.add_vertex(rev_complement(kmer2))
            self.add_edge(rev_complement(kmer2), rev_complement(kmer1))


    def get_contig_wrapper(self, start_node, visit):

        def is_self_ring(v):
            return v.seq in v.out_edge

        def select_valid_path(v):
            """
            select a valid path to traverse
            """
            if is_self_ring(v):
                return None
            cmpt_edge = list(v.out_edge.values())
            # un-ambiguous path
            if len(cmpt_edge) == 1:
                return cmpt_edge[0]
            
            elif len(cmpt_edge) >= 2:
                cmpt_edge = sorted(cmpt_edge, key=lambda x:x.cov)
                if cmpt_edge[-1].cov > cmpt_edge[-2].cov:
                    return cmpt_edge[-1]
                       

        node_stack = collections.deque()
        node_stack.append(start_node)
        contigs = []
        a = 0
        # stack stores nodes whose degree != 2
        while len(node_stack):
            a += 1
            b = 0
            node = node_stack.pop()
            visit.add(node.seq)
            flg = True
            #print(len(self.v), len(visit), a)
            for edge in node.out_edge.values():
                # travel until degree > 2 or degree = 1
                prefix = ''
                path = []
                v = self.v[edge.seq2]
                while True:
                    b += 1
                    path.append(v.seq[-1])
                    if v.seq in visit:
                        break
                    visit.add(v.seq)
                    sel_edge = select_valid_path(v)
                    if sel_edge is None:
                        if len(v.out_edge) >= 2:
                            node_stack.append(v)
                        break
                    else:
                        v = self.v[sel_edge.seq2]
                contigs.append(prefix + ''.join(path))
                flg = False
        return contigs           


    def get_contig_graph(self):
        """
        dfs and store any unambiguous contig
        """
        visit = set()
        all_path = []
        c = 0
        """
        while True:
            c += 1
            print(c,len(self.v),len(visit))
            if c > 2000:
                break
            start_node = None
            for k in self.v:
                if k not in visit and len(self.v[k].out_edge)!=2:
                    start_node = self.v[k]
                    break
            if not start_node:
                break
            paths = self.get_contig_wrapper(start_node, visit)
            all_path.extend(paths)
        """
        for v in self.v:
            if v not in visit and len(self.v[v].out_edge) != 2:
                start_node = self.v[v]
                paths = self.get_contig_wrapper(start_node, visit)
                all_path.extend(paths)
                c += 1
                print(c,len(self.v),len(visit))
        return all_path

def write_fa(f, lines):
    for i,line in enumerate(lines):
            if len(line.strip()) > 50:
                f.write('>{} length {} xxx\n'.format(i*2+1, len(line.strip())))
                f.write(line.strip()+'\n')


if __name__ == '__main__':
    g = DBG(k=51)
    reads = read_fasta('data/data3/short_1.fasta')
    reads += read_fasta('data/data3/short_2.fasta')
    #reads = ['TTCTACTATCGCTGTGGGATGGATCATAAA',
           #  'TTCTACTATCGCTGTGGGATGGATCATCCC',
    #         'AAATTCCCCCCCCCCCCCC']
    for i,read in enumerate(reads):
        g.add_read(read)
        print(i)
        #if i == 100:
        #    break
    contigs = g.get_contig_graph()
    #contigs = g.dfs_graph()
    f = open('result/data3_contig.txt','w')
    write_fa(f, contigs)