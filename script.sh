#!/bin/sh

contig_path="result/data2_contig17.txt"
data_dir="data/data2"

cd ./external/dbg2olc/
./DBG2OLC k 13 AdaptiveTh 0.002 KmerCovTh 2 MinOverlap 5 RemoveChimera 0 Contigs ../../${contig_path} f ../../${data_dir}/long.fasta 
cd ../../
cp ./external/dbg2olc/DBG2OLC_Consensus_info.txt ./external/sparc
cp ./external/dbg2olc/backbone_raw.fasta ./external/sparc
cat ${contig_path} ${data_dir}/long.fasta > ./external/sparc/all_ctg.fasta
cd ./external/sparc
./split_and_run_sparc.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt all_ctg.fasta ./consensus_dir
cd ../../
cat ./external/sparc/consensus_dir/final_assembly.fasta ${contig_path} > result/full2_contig_2.txt
