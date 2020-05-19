```python
######################################################################################################
def rndTmpFileName(tmpDirName='/tmp'):
  return tmpDirName+'/'+''.join(random.choices(string.digits,k=20))

######################################################################################################
def buildGas(gffFp=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  genes=HTSeq.GenomicArrayOfSets("auto",stranded=False)
  for feature in HTSeq.GFF_Reader(str(gffFp)):
    if feature.type in ('gene'):
      genes[feature.iv]+=(feature.attr['ID'],feature.iv.strand)
  return genes

######################################################################################################
def starFunc(args):
  if args[-1]: wlog(functionCallInfo(frame=inspect.currentframe()))
  out=(lambda x:x[0](*x[1:]))(args)
  wlog(out)
  return out

######################################################################################################
def parallelGetAlns(inputBamFpsDict=None,featsGas=None,outputGffAlnsFpDict=None,cpu=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  if not all(fp.is_file() and fp.stat().st_size>0 for fp in outputGffAlnsFpDict.values()):
    argsList=\
      (lambda ll:random.sample(ll, k=len(ll)))\
      (list(itertools.chain.from_iterable(\
        list((getHybridSortedPairRanks,k,v,refId,workDirPaths['tmp']) for refId in pysam.AlignmentFile(v).references) \
        for k,v in inputBamFpsDict.items())))
    wlog('parallelGetFaAlns: multiprocessing getHybridSortedPairRanks')
    with multiprocessing.Pool(processes=cpu) as pool:
      compressedHybridSortedPairRanksDict=dict(\
        (k,list(e[1:] for e in g)) \
        for k,g in itertools.groupby(sorted(pool.imap_unordered(starFunc,argsList),key=lambda x:(x[0],x[1])),key=lambda x:x[0]))
      pool.close() ; pool.join()
    wlog('compressedHybridSortedPairRanksDict created')
    argsList=list((getGenomicAlns,k,inputBamFpsDict[k],v,featsGas,outputGffAlnsFpDict[k],True) for k,v in compressedHybridSortedPairRanksDict.items())
    wlog('parallelGetFaAlns: multiprocessing getGenomicAlns')
    with multiprocessing.Pool(processes=cpu) as pool:
      results=dict((tag,outputGffAlnsFp) for tag,outputGffAlnsFp in pool.imap_unordered(starFunc,argsList))
      pool.close() ; pool.join()
    wlog('outputGffAlnsFpDict completed : '+str(all(outputGffAlnsFpDict[tag]==outputGffAlnsFp for tag,outputGffAlnsFp in results.items())))
  return None if not all(fp.is_file() and fp.stat().st_size>0 for fp in outputGffAlnsFpDict.values()) else outputGffAlnsFpDict

######################################################################################################
def getHybridSortedPairRanks(tag=None,inputBamFp=None,refId=None,tmpDir=None):
  pairRanksFileName=rndTmpFileName(tmpDirName=str(tmpDir))
  rankedSamIterator=enumerate(pysam.AlignmentFile(inputBamFp).fetch(refId))
  with open(pairRanksFileName,'wb') as fh:
    (lambda compressor,ll:(fh.write(compressor.compress(struct.pack(str(len(ll))+'I',*ll))+compressor.flush()),fh.flush(),fh.seek(0)))\
    (zlib.compressobj(level=1,wbits=-15),list(itertools.chain(*sorted(\
      (lambda tt:(tt[0],tt[1]-tt[0]) if len(tt)>1 else (tt[0],0))(list(zip(*g))[2]) \
      for k,g in itertools.groupby(sorted(iter((aln.query_name,aln.reference_start,e) for e,aln in rankedSamIterator),key=lambda x:(x[:2])),key=lambda x:x[0])))))
  return tag,refId,pathlib.Path(pairRanksFileName)

######################################################################################################
def bamHybridSort(tag=None,pairRanksFp=None,rankedSamIterator=None):
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
  with pairRanksFp.open('rb') as fh:
    decompressor=zlib.decompressobj(wbits=-15) ; decompData=b'' ; rank=-1
    chunk=fh.read(2048)
    while chunk:
      decompData+=decompressor.decompress(chunk)
      decompData,rankedSamIterator,rank,alnsPairsList=processDecompData(decompData=decompData,rankedSamIterator=rankedSamIterator,rank=rank)
      for alnsPair in alnsPairsList: yield alnsPair
      chunk=fh.read(2048)
    decompData+=decompressor.flush()
    decompData,rankedSamIterator,rank,alnsPairsList=processDecompData(decompData=decompData,rankedSamIterator=rankedSamIterator,rank=rank)
    for alnsPair in alnsPairsList: yield alnsPair
    assert decompData==b''
  pairRanksFp.unlink()

######################################################################################################
def getGenomicAlns(tag=None,inputBamFp=None,pairRanksFpsList=None,featsGas=None,outputGffAlnsFp=None,log=True):
  if log: wlog(functionCallInfo(frame=inspect.currentframe()))
  ######################################################################################################
  def getNtPosList(aln=None):
    return \
      (lambda ll:list(zip(*((pos,nt) for pos,nt in functools.reduce(lambda x,y:x+[y] if y[0] else x[:-1]+[(x[-1][0],x[-1][1]+y[1])],ll[1:],[ll[0]])))))\
      (list((p[1],aln.query_alignment_sequence[p[0]] if p[0]!=None else 'N') for p in aln.get_aligned_pairs()))
  ######################################################################################################
  def getGffDict(readsBlock=None,refId=None,strand=None):
    ll1,ll2,ll3=list(zip(*(\
      (lambda tt,T:(pos,tt[0] if (tt[1]/T)>0.5 else 'X',tt[1] if not(tt[0] in ('N','X')) else 0))\
      (counter.most_common(1)[0],sum(counter.values())) for pos,counter in readsBlock)))
    ll4=list(itertools.chain(*([(1,ntStr[0] if ntStr[0] in ('N','X') else 'M')]+list((1,'I') for nt in ntStr[1:]) for ntStr in ll2))) ; revertSign=(strand=='+')*2-1
    gffDict=collections.OrderedDict(\
      [('seqid',refId),('source','pisumFrostRnaSeq'),('type','match'),('start',ll1[0]),('end',ll1[-1]),('score','.'),('strand',strand),('phase','.')]+\
      [('attributes',collections.OrderedDict(\
        [('ID',''),('seq',revComp(inputSeq=''.join(ll2),rev=(strand=='-'),comp=(strand=='-'))),('coverage',round(sum(ll3)/len(list(filter(lambda x:x,ll3))),2)),('ovlFeatIds','')]+\
        [('Gap',''.join(''.join(map(str,o)) for o in functools.reduce(lambda x,y:x[:-1]+[(x[-1][0]+1,x[-1][1])] if y[1]==x[-1][1] else x+[y],ll4[1:],[ll4[0]])[::revertSign]))]))])
    return gffDict
  ######################################################################################################
  def fileWriteGffAndFaAlns(gffDict=None,featsGas=None,outputGffAlnsFh=None,outputFaAlnsFh=None,cnt=None):
    cnt+=1
    if refId in featsGas.chrom_vectors:
      gffDictIv=HTSeq.GenomicInterval(chrom=refId,start=gffDict['start'],end=gffDict['end'],strand=gffDict['strand'])
      ovlFeatIdsSet='|'.join(sorted(set(featId[0] for iv,featIdsSet in featsGas[gffDictIv].steps() for featId in featIdsSet if featId and featId[1]==gffDict['strand'])))
    else: ovlFeatIdsSet=''
    gffDict['attributes']['ID']='match'+'%06d'%cnt+'_'+refId
    gffDict['attributes']['ovlFeatIds']=ovlFeatIdsSet or 'None'
    attributes=gffDict.pop('attributes')
    tmpOutputFaAlnsFh.write(fastaFormat(seqId=attributes['ID'],inputSeq=attributes.pop('seq'))+'\n')
    outputGffAlnsFh.write('\t'.join(map(str,gffDict.values()))+'\t'+';'.join(k+':'+str(v) for k,v in attributes.items())+'\n')
    outputGffAlnsFh.flush() ; tmpOutputFaAlnsFh.flush()
    return cnt
  if not(outputGffAlnsFp.is_file() and outputGffAlnsFp.stat().st_size>0):
    refIdsList=pysam.AlignmentFile(inputBamFp).references ; flagToStrandDict={0:'-',73:'-',99:'-',147:'-',153:'-',16:'+',83:'+',89:'+',137:'+',163:'+'} ; strands=('-','+')
    with outputGffAlnsFp.open('wt') as outputGffAlnsFh, tempfile.NamedTemporaryFile('w+t') as tmpOutputFaAlnsFh:
      for refId,pairRanksFp in pairRanksFpsList:
        wlog('getting matches on '+refId+' for '+tag)
        dictOfNtPosLists=dict((strand,list()) for strand in strands) ; gffDictsDeque=collections.deque() ; cnt=0
        for alnsPair in bamHybridSort(tag=tag,pairRanksFp=pairRanksFp,rankedSamIterator=enumerate(pysam.AlignmentFile(inputBamFp))):
          strand=flagToStrandDict[alnsPair[0].flag]
          ntPosList=list(list(itertools.chain(*ll)) for ll in zip(*(getNtPosList(aln=aln) for aln in alnsPair)))
          if dictOfNtPosLists[strand]==[]:
            dictOfNtPosLists[strand]=[ntPosList[0][-1],ntPosList]
          elif ntPosList[0][0]<=dictOfNtPosLists[strand][0]:
            for i in (0,1): dictOfNtPosLists[strand][1][i].extend(ntPosList[i])
            dictOfNtPosLists[strand][0]=max(dictOfNtPosLists[strand][0],ntPosList[0][-1])
          else:
            readsBlock=(lambda dd1:(list(dd1[p].update({n:1}) for p,n in zip(*dictOfNtPosLists[strand][1])),dd1)[1])(collections.defaultdict(collections.Counter))
            readsBlock=list((pos,readsBlock[pos] if pos in readsBlock.keys() else collections.Counter({'X':1})) for pos in range(min(readsBlock.keys()),max(readsBlock.keys())+1))
            dictOfNtPosLists[strand]=[ntPosList[0][-1],ntPosList]
            gffDictsDeque.append(getGffDict(readsBlock=readsBlock,refId=refId,strand=strand))
            gffDictsDeque=collections.deque(sorted(gffDictsDeque,key=lambda dd:dd['start']))            
            while len(set(dd['strand'] for dd in gffDictsDeque))>1:
              cnt=fileWriteGffAndFaAlns(gffDict=gffDictsDeque.popleft(),featsGas=featsGas,outputGffAlnsFh=outputGffAlnsFh,outputFaAlnsFh=tmpOutputFaAlnsFh,cnt=cnt)
        for strand in strands:
          if dictOfNtPosLists[strand]:
            readsBlock=(lambda dd1:(list(dd1[p].update({n:1}) for p,n in zip(*dictOfNtPosLists[strand][1])),dd1)[1])(collections.defaultdict(collections.Counter))
            readsBlock=list((pos,readsBlock[pos] if pos in readsBlock.keys() else collections.Counter({'X':1})) for pos in range(min(readsBlock.keys()),max(readsBlock.keys())+1))
            gffDictsDeque.append(getGffDict(readsBlock=readsBlock,refId=refId,strand=strand))
        gffDictsDeque=collections.deque(sorted(gffDictsDeque,key=lambda dd:dd['start']))
        for gffDict in gffDictsDeque:
          cnt=fileWriteGffAndFaAlns(gffDict=gffDict,featsGas=featsGas,outputGffAlnsFh=outputGffAlnsFh,outputFaAlnsFh=tmpOutputFaAlnsFh,cnt=cnt)
      tmpOutputFaAlnsFh.seek(0) ; outputGffAlnsFh.write('##FASTA'+'\n'+tmpOutputFaAlnsFh.read().strip())
  return None if not(outputGffAlnsFp.is_file() and outputGffAlnsFp.stat().st_size>0) else tag,outputGffAlnsFp

  ```
