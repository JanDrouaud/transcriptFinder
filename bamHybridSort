```ruby

import collections, io, itertools, pathlib, pysam, struct, zlib

######################################################################################################
def bamHybridSort(pysamIterator=None):
  ######################################################################################################
  def getPairRanksStream(rankedSamIterator=None):
    # ~ # sort alns by name then position, appending rank from the initial file to each aln
    ll1=collections.deque(\
      sorted(iter((aln.query_name,aln.reference_start,e) for e,aln in rankedSamIterator),key=lambda x:(x[:2])))
    # If no aln in iterator, return a void stream
    if not ll1: return io.BytesIO()
    # group into pairs
    ll2=[[ll1.popleft()]]
    while ll1:
      x=ll1.popleft()
      if x[0]==ll2[-1][0][0]:ll2[-1]+=[x]
      else: ll2+=[[x]]
    # For each pair of alns, extract the rank of the first one and compute the offset to the rank of the second one
    # (0 if not a pair). Chain all pairs of values then write compressed data to a stream
    return \
      (lambda stream,compressor,ll:(stream.write(compressor.compress(struct.pack(str(len(ll))+'I',*ll))+compressor.flush()),stream.flush(),stream.seek(0),stream)[3])\
      (io.BytesIO(),zlib.compressobj(level=1,wbits=-15),list(itertools.chain(*sorted(iter((lambda p2:[p2[0],len(p2)>1 and p2[1]-p2[0]])(list(a[2] for a in p1)) for p1 in ll2)))))
  ######################################################################################################
  def processDecompData(decompData=None,rankedSamIterator=None,rank=None):
    alnsPairsList=list()
    # each alnRanksPair require 8 bytes (two 4 bytes integer)
    alnsPairsNumber=len(decompData)//8
    # process the pairs of alns, yielding each aln
    for alnRanksPair in list(itertools.zip_longest(*[iter(struct.unpack(str(alnsPairsNumber*2)+'I',decompData[:alnsPairsNumber*8]))]*2)):
      while rank<alnRanksPair[0]: rank,aln1=next(rankedSamIterator)
      if not(alnRanksPair[1]): alnsPairsList.append((aln1,))
      else:
        rankedSamIterator,it=list(itertools.tee(rankedSamIterator,2))
        for i in range(alnRanksPair[1]): f,aln2=next(it)
        alnsPairsList.append((aln1,aln2))
    return decompData[alnsPairsNumber*8:],rankedSamIterator,rank,alnsPairsList
  rankedSamIterator1,rankedSamIterator2=list(map(enumerate,itertools.tee(pysamIterator,2)))
  # ~ (lambda it1,it2:(enumerate(it1),it2))(*itertools.tee(samIterator,2))
  pairRanksStream=getPairRanksStream(rankedSamIterator=rankedSamIterator1) ; rank=-1 ; decompressor=zlib.decompressobj(wbits=-15) ; decompData=b''
  chunk=pairRanksStream.read(2048)
  while chunk:
    decompData+=decompressor.decompress(chunk)
    decompData,rankedSamIterator2,rank,alnsPairsList=processDecompData(decompData=decompData,rankedSamIterator=rankedSamIterator2,rank=rank)
    for alnsPair in alnsPairsList: yield alnsPair
    chunk=pairRanksStream.read(2048)
  decompData+=decompressor.flush()
  decompData,rankedSamIterator2,rank,alnsPairsList=processDecompData(decompData=decompData,rankedSamIterator=rankedSamIterator2,rank=rank)
  for alnsPair in alnsPairsList: yield alnsPair
  assert decompData==b''
  pairRanksStream.close()
  
  ```
