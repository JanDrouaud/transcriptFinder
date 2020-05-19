# A set of memory economical functions for fast parallel assembly of genome mapped RNA-seq reads

Well, say we have to cope with a very large genome, human for example. Say we have a set of bam files holding alignments of stranded RNA-seq reads to the genome sequence. We want to dig out transcriptomic data, that is the **sequence** of the transcripts which are actually present in the sample tissue/organ/organism of interest.

Usually, alignments in bam files are sorted either by reference (i.e. pseudo-chromosome/scaffold/contig) ID, then coordinate, or by read name. For the purpose of grouping together all alignments that map "overlappingly" to the genome sequence, we have first to perform a kind of **hybrid sorting**: alignments are sorted by reference id / coordinate only if their come single-end read or first member of paired-ends reads. The second member of paired-ends reads, that map farther than the first one, is placed - in that hybrid sorting method - immediately next to its kin.

Blocks of alignments - from reads that originate from the **same strand** - are then all aligned together, simply considering the genomic coordinates Then a consensus sequence and associated CIGAR line are computed. If the mapping (stranded) region overlaps with some known genomic features, their ID is also provided. Yeahh, parsing all that data in clean gff3 compliant files, we can visualize it with most genome browsers.

Presently, only the consensus sequence is extracted, choosing at every genomic position the most common nucleotide. I may consider in some future to build every possible transcript, using a graph-based method. Is it really a good idea ?
