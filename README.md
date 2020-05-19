# A set of memory economical functions for fast parallel assembly of genome mapped RNA-seq reads

Well, say we have to cope with a very large genome, human for example.

Given RNA-seq data in bam files holding both single-end and paired-ends alignments, sorted by chromosome/scaffold then coordinate, we want to iterate over alignment pairs, sorted by first alignment coordinate.
The bamHybridSort function has been designed for that purpose. It uses a compressed data stream in order to reduce the memory workload, which can get huge when processing large chromosomes.
