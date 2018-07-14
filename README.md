# Genome Assembly Implementation

This is a course project for Algorithm Design and Analysis in Fudan University.

The pipeline follows paper "Hybrid error correction and de novo assembly of single-molecule sequencing reads".

Outline:
- Use DBG(de bruijn graph) to analyze reads
- Heuristic error correction
- Take reliable contigs as input for OLC(over lap consensus) solver(external implementation [here](https://github.com/yechengxi/DBG2OLC))

**This repository holds personal implementation. The code is helpful to algorithm learners, but for advanced usage of biomedical research, please refer to https://github.com/yechengxi/DBG2OLC**
