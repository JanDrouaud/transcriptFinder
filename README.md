# A set of memory economical functions for fast parallel assembly of genome mapped RNA-seq reads

Well, say we have to cope with a very large genome, human for example. Say we have a set of bam files holding alignments of stranded RNA-seq reads to the genome sequence. We want to dig out transcriptomic data, that is the **sequence** of the transcripts which are actually present in the sample tissue/organ/organism of interest.

Usually, alignments in bam files are sorted either by reference (i.e. pseudo-chromosome/scaffold/contig) id, then coordinate, or by read name. For the purpose of grouping together all alignments that map *overlappingly* to the genome sequence, we have first to perform a kind of *hybrid* sorting: alignments are sorted by reference id / coordinate only if their come single-end read or first member of paired-ends reads. The second member of paired-ends reads, that map farther than the first one, is placed - in that hybrid sorting method - immediately next to its kin.


The set of functions in 
Given RNA-seq data in bam files holding both single-end and paired-ends alignments, sorted by chromosome/scaffold then coordinate, we will iterate over alignment pairs, sorted by first alignment coordinate. All reads
The bamHybridSort function has been designed for that purpose. It uses a compressed data stream in order to reduce the memory workload, which can get huge when processing large chromosomes.
