- [X] Commit 5ccc063 of transcriptFinder.py is not memory economical: there is a huge peak of memory use ending in a 'Memory Error' crash. This originate from a transcriptomic peculiarity : an assembly of several 10^5 alignments, around 20 Kb wide and 35000 mean depth-of-coverage. To fix urgently.
- [ ] Previous bug has been fixed in commit 7a8258cc. The getGenomicAlns function has been slightly modified, so as to use less memory. Consequently, the running time increases slightly (10 % ?). To optimize.