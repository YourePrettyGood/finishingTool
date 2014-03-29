
import os
from itertools import groupby
from collections import deque
from operator import itemgetter
import sys
import time

##### House keeping files
def writeToFile_Double1(folderName, fileName1, fileName2, option = "contig"):

    f2 = open(folderName + fileName2, 'w')
    fOriginal = open(folderName + fileName1, 'r')
    
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) >0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead+ tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)
    
    print "len(readSet)", len(readSet)
    
    fOriginal.close()
    
    if option == "contig":
        header = ">Contig"
    else:
        header = ">Read"
    for eachcontig, dum in zip(readSet, range(len(readSet))):
        f2.write(header+ str(dum)+"_p\n")
        f2.write(eachcontig+'\n') 
        f2.write(header+ str(dum)+"_d\n")
        f2.write(reverseComplement(eachcontig)+'\n')
        
    f2.close()
                   
   
def writeToFile_Double2(folderName, fileName1, fileName2):
    # for reads
    
    f = open(folderName + fileName1, 'r')
    f2 = open(folderName + fileName2, 'w')
    
    tmp1 = f.readline().rstrip()
    if len(tmp1) > 0:
        tmp2 = f.readline().rstrip()
        
    while len(tmp1) > 0 and len(tmp2) >0:
        f2.write(tmp1+ "_p")
        f2.write('\n')
        f2.write(tmp2)
        f2.write('\n')
        f2.write(tmp1 + "_d")
        f2.write('\n')
        f2.write(reverseComplement(tmp2))
        f2.write('\n')
        
        tmp1 = f.readline().rstrip()
        
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
    
    f2.close()
    f.close()

def writeToFile(f2, runningIndex,seq):
    f2.write(">Seg_"+ str(runningIndex))
    f2.write('\n')
    f2.write(seq)
    f2.write('\n')
   

def findContigLength(folderName, option):
    
    if option == "contigs":
        print "\n\nfindContigLength(folderName)\n"
        contigLength = {}
        f = open(folderName + "contigs_Double.fasta", 'r')
        
        
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) >0:
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()
        
        
    else:
        print "\n\nfindContigLength(folderName)\n"
        contigLength = {}
        f = open(folderName + "improved_Double.fasta", 'r')
        
        
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) >0:
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()
        
    
    
        f = open(folderName + "relatedReads_Double.fasta", 'r')
        
    
        tmp1 = f.readline().rstrip()
        if len(tmp1) > 0:
            tmp2 = f.readline().rstrip()
            
        while len(tmp1) > 0 and len(tmp2) >0:
    
            contigLength[tmp1[1:]] = len(tmp2)
            
            tmp1 = f.readline().rstrip()
            
            if len(tmp1) > 0:
                tmp2 = f.readline().rstrip()
                
        f.close()

    return contigLength

def reverseComplement(myStr):
    myNewStr = myStr[::-1]
    myNewStr2 = ""
    for i in range(len(myNewStr)):
        if myNewStr[i] == 'A' or myNewStr[i] == 'a':
            myNewStr2 += 'T'
            
        elif  myNewStr[i] == 'T' or myNewStr[i] == 't':
            myNewStr2 +='A'
            
        elif  myNewStr[i] == 'C' or myNewStr[i] == 'c':
            myNewStr2 += 'G'
            
        elif  myNewStr[i] == 'G' or myNewStr[i] == 'g':
            myNewStr2 += 'C'
            
    return myNewStr2

##### New method of using string graph to approximate DB graph

    
def formRelatedReadsFile(folderName,mummerLink):    
    # Find associated read and extract into a file associatedReads.fasta
    # Input: contigs.fasta, cleaned_Reads.fasta 
    # Output: relatedReads.fasta

    ### Extract heads of the contigs
    print ">formRelatedReadsFile"
    f = open(folderName + "improved.fasta", 'r')
    f2 = open(folderName + "improvedTrunc.fasta", 'w')
    temp = f.readline()
    tempContig = ""
    thres = 400
    runningIndex = 0
    endThres = 10 
    if True:
        while len(temp) > 0:
            if temp[-1] == '\n':
                temp = temp[0:-1]
            
            
            if temp[0] == '>':
    
                if len(tempContig) > 0:
                    writeToFile(f2, runningIndex,tempContig[0:thres])
                    runningIndex = runningIndex +1
                    
                    writeToFile(f2, runningIndex,tempContig[-thres:])
                    runningIndex = runningIndex +1 
                    
                                    
                    writeToFile(f2, runningIndex,reverseComplement(tempContig[0:thres]))
                    runningIndex = runningIndex +1
                    
                    writeToFile(f2, runningIndex,reverseComplement(tempContig[-thres:]))
                    runningIndex = runningIndex +1
                    
                    tempContig = ""
            else:
                tempContig = tempContig + temp
            
            temp = f.readline()
    
        writeToFile(f2, runningIndex,tempContig[0:thres])
        runningIndex = runningIndex +1
        
        writeToFile(f2, runningIndex,tempContig[-thres:])
        runningIndex = runningIndex +1
        
                        
        writeToFile(f2, runningIndex,reverseComplement(tempContig[0:thres]))
        runningIndex = runningIndex +1
        
        writeToFile(f2, runningIndex,reverseComplement(tempContig[-thres:]))
        runningIndex = runningIndex +1
        
        
        f2.close()
        f.close()
        
        ### Write double stranded reads
        writeToFile_Double1(folderName, "improved.fasta", "improved_Double.fasta", "contig")
        writeToFile_Double1(folderName, "raw_reads.fasta", "raw_reads_Double.fasta","read")
    
        
        ### Apply MUMMER on them using cleanedReads against them
        command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"out "+ folderName+ "improvedTrunc.fasta "+ folderName+ "raw_reads.fasta"
        os.system(command)
        
        command  = mummerLink +"show-coords -r "+folderName+"out.delta > "+folderName+"fromMum"
        os.system(command)
    
    f = open(folderName + "fromMum", 'r')
    assoiatedReadIndex = []
    nameList = []
    for i in range(6):
        tmp = f.readline()
    
    while len(tmp) > 0:
        infoArr = tmp.split('|')
        myArr = infoArr[-1].split('\t')
        rdGpArr = infoArr[-1].split('\t')
        contigName = rdGpArr[0].rstrip().lstrip()
        readName = rdGpArr[1].rstrip().lstrip()
            
        
        
        endSegArr = infoArr[0].split(" ")
        pos = []
        for eachitem in endSegArr:
            if len(eachitem) > 0:
                pos.append(int(eachitem))
                
        startPos = pos[0]
        endPos = pos[1]
        if startPos < endThres and endPos > thres-endThres:
            assoiatedReadIndex.append(myArr[1])
            nameList.append([int(contigName.split('_')[1]), readName])
        tmp  = f.readline()
    
    f.close()
    
    
    nameList.sort()

    assoiatedReadIndex.sort()
    
    print "assoiatedReadIndex", assoiatedReadIndex
    
    ckIndex = 0
    f = open(folderName+"associatedNames.txt", 'w')
    oneItem = 0
    keyFound = []
    for key, items in groupby(assoiatedReadIndex):
        
        countItem = 0
        for eachitem in items:
            countItem += 1
            

        if countItem == 1:
            
            oneItem += 1
        else:
            key = key.rstrip()
            if not key in keyFound:
                f.write(key+'\n')
                keyFound.append(key)

        ckIndex += 1
    
    print "ckIndex,oneItem: ",ckIndex, oneItem
    f.close()

    command = "perl -ne 'if(/^>(\S+)/){$c=$i{$1}}$c?print:chomp;$i{$_}=1 if @ARGV' "+folderName+"associatedNames.txt "+folderName+"raw_reads.fasta > "+folderName+"relatedReads.fasta"
    os.system(command)
    
    writeToFile_Double1(folderName, "relatedReads.fasta", "relatedReads_Double.fasta","read")

    
    
    
def extractEdgeSet(folderName, mummerLink):
    # Tasks: reconstruct the string  graph
    
    # Input : relatedReads_Double.fasta, conig_Double.fasta
    # Intermediate files: fromMum_overlap , fromMum_overlap
    # Output: connectivity of eachNode: InList, OutList [critical]
    #         connectivity of eachNode: arrow representatiowith size [optional]
    
    
    ### Perform MUMMER alignment
    print ">Extract Edge set"
    lengthDic = findContigLength(folderName, "improved")
    numberOfContig = 0
    f = open(folderName+ "improved_Double.fasta",'r')
    tmp = f.readline()
    tmp = tmp.rstrip()
    while (len(tmp)>0):
        numberOfContig += 1
        tmp = f.readline()
        tmp = tmp.rstrip()
    
    numberOfContig = numberOfContig /2
    print "numberOfContig", numberOfContig

    f.close()
    K = 200
    
    thres = 7
    # Nodes are contigs, 
    nodes = [i for i in range(numberOfContig)]
    dataSet = []
    
    ### Apply MUMMER on them using cleanedReads against them
    command =mummerLink +"nucmer --maxmatch --simplify -p "+folderName+"outRefine "+ folderName+ "improved_Double.fasta "+ folderName+ "relatedReads_Double.fasta"
    os.system(command)
    
    command  = mummerLink +"show-coords -r "+folderName+"outRefine.delta > "+folderName+"fromMumRefine"
    os.system(command)
    

    
    f = open(folderName + "fromMumRefine", 'r')
    for i in range(6):
        tmp = f.readline()

    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr =  info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        matchLenArr = info[2].split()
       
    
        matchLen = int(matchLenArr[0])    
        contigStart, contigEnd =  int( firstArr[0]), int( firstArr[1])
        readStart, readEnd =  int(filterArr[0]) , int(filterArr[1])
        
        contigName = rdGpArr[0].rstrip().lstrip()
        readName = rdGpArr[1].rstrip().lstrip()
            
    
        if readStart < readEnd and matchLen> K and min(contigStart,readStart)  < thres and min(lengthDic[contigName]- contigEnd,  lengthDic[readName] - readEnd) < thres:
            conditionForMatch = True
        else:
            conditionForMatch = False

        if conditionForMatch :
            if contigStart < thres:
                dataSet.append((readName, contigName, 'L',matchLen))
            
            if lengthDic[contigName]- contigEnd < thres :
                dataSet.append((readName, contigName, 'R',matchLen))
        
        tmp = f.readline()
        
    f.close()  
    dataSet.sort()
    
    print "dataSet:", dataSet[0]

    matchPair = []
    for key, items in groupby(dataSet, itemgetter(0)):
        left = []
        right = []
        
        for subitem in items:

            myArr =subitem[1].split('_')
            orientation = myArr[1]

            if orientation == 'p' :
                contigNum = int(myArr[0][6:])*2
            else:
                contigNum = int(myArr[0][6:])*2 +1
            
            print subitem[3]
            if subitem[2] == 'L':
                left.append([contigNum, subitem[3]])
            else:
                right.append([contigNum, subitem[3]])
                
        print left, right
        #if len(set(left).intersection(set(right))) == 0 and len(set(left))== len(left) and len(set(right)) == len(right):
        #if  len(set(left))== len(left) and len(set(right)) == len(right
            
        for eachleft in left:
            for eachright in right:
                leftIndex , rightIndex = eachleft[0], eachright[0]
                leftLen, rightLen = eachleft[1], eachright[1]
                
                if leftIndex != rightIndex:
                    matchPair.append([rightIndex, leftIndex, min(leftLen,rightLen)])

    matchPair.sort()
        
    keyFound = []
    bestMatchPair = []
    for key, items in groupby(matchPair, itemgetter(0,1)):
        maxvalue = -1
        for eachitem in items:
            if eachitem[2] > maxvalue:
                maxvalue = eachitem[2]
        bestMatchPair.append([key[0], key[1], maxvalue])
    
    bestMatchPair.sort( key=lambda tup: tup[2], reverse=True)
    
    print bestMatchPair


    leftConnect = [-1 for i in range(numberOfContig)]
    rightConnect = [-1 for i in range(numberOfContig)]
    
    for eachitem in bestMatchPair:

        prefixContig = eachitem[0]
        suffixContig = eachitem[1]
        
        if leftConnect[suffixContig] == -1 and rightConnect[prefixContig] == -1:
            leftConnect[suffixContig] = prefixContig 
            rightConnect[prefixContig] = suffixContig
    
    startList= []
    print leftConnect
    for i in range(len(leftConnect)):
        if leftConnect[i] == -1:
            startList.append(i)
    
    contigList = []
    print startList
    for eachitem in startList:
        tmp = eachitem
        myList = [tmp]
        while rightConnect[tmp] != -1:
            tmp = rightConnect[tmp]
            if tmp != -1:
                myList.append(tmp)
        contigList.append(myList)
    
    print contigList

            
    i = 0
    fOriginal = open(folderName + "improved.fasta", 'r')
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) >0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead+ tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)  
    
    seqToPrint = []
    contigUsed = [False for i in range(numberOfContig/2)]
    fAllInOne = open(folderName+"allInOne2.fasta",'w' )
    finalList = []
    for eachContig, i in zip(contigList, range(len(contigList))):
        tmpList = []
        for eachitem in eachContig:

            readNum = eachitem/2
            if contigUsed[readNum] ==False:
                seqToPrint.append(eachitem)
                tmpList.append(eachitem)
                contigUsed[readNum] = True
        if len(tmpList) > 0:
            finalList.append(tmpList)
        
    
    fImproved = open(folderName +"improved2.fasta", 'w')
    for eachcontig, dummyIndex in zip(finalList, range(len(finalList))):
        fImproved.write(">Segkk"+str(dummyIndex)+'\n')
        for eachseg in eachcontig:
            readNum = eachseg/2
            orientation = eachseg%2

            if orientation == 0:
                fImproved.write(readSet[readNum])
            else:
                fImproved.write(reverseComplement(readSet[readNum]))
        fImproved.write('\n')
        
    fImproved.close()
            


    for eachContigIndex,dum in zip(seqToPrint,range(len(seqToPrint))):

        readNum = eachContigIndex/2
        orientation = eachContigIndex%2
        
        fAllInOne.write(">Seg_"+str(dum)+'\n')
        if orientation == 0:
            fAllInOne.write(readSet[readNum])
            fAllInOne.write('\n')
        else:
            fAllInOne.write(reverseComplement(readSet[readNum]))
            fAllInOne.write('\n')
                    
    fAllInOne.close()
    


    for eachcontig, i in zip(finalList, range(len(finalList))):
        fout = open(folderName + "improvedContig2_"+str(i)+".fasta", 'w')
        seqToPrint = []

        for eachitem in eachcontig:

            readNum = eachitem/2
            seqToPrint.append(eachitem)

        print "ImprovedContig ",i 
        for eachhaha in seqToPrint:
            print len(readSet[eachhaha/2])
            
        for eachContigIndex,dum in zip(seqToPrint,range(len(seqToPrint))):

            readNum = eachContigIndex/2
            orientation = eachContigIndex%2
            
            fout.write(">Seg_"+str(dum)+'\n')
            if orientation == 0:
                fout.write(readSet[readNum])
                fout.write('\n')
            else:
                fout.write(reverseComplement(readSet[readNum]))
                fout.write('\n')


        i += 1
        fout.close()




    
def newGraphPipeLine(folderName, mummerLink):
    print "newGraphPipeLine"
    formRelatedReadsFile(folderName,mummerLink)
    extractEdgeSet(folderName, mummerLink)
    
    
def greedyAlg(mummerLink, folderName):
    
    print "Direct greedy"
    
    thres = 7
    writeToFile_Double1(folderName, "contigs.fasta", "contigs_Double.fasta", "contig")
    
    command =mummerLink +"nucmer --maxmatch --nosimplify -p "+folderName+"greedy "+ folderName+ "contigs_Double.fasta "+ folderName+ "contigs_Double.fasta"
    os.system(command)
    
    command  = mummerLink +"show-coords -r "+folderName+"greedy.delta > "+folderName+"fromMumGreedy"
    os.system(command)
        
    lengthDic = findContigLength(folderName, "contigs")
    

    f = open(folderName + "fromMumGreedy", 'r')
    for i in range(6):
        tmp = f.readline()
        
    dataSet = []
    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr =  info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        matchLenArr = info[2].split()
       
    
        matchLen = int(matchLenArr[0])    
        contigStart, contigEnd =  int( firstArr[0]), int( firstArr[1])
        contig2Start, contig2End =  int(filterArr[0]) , int(filterArr[1])
        
        contigName = rdGpArr[0].rstrip().lstrip()
        contig2Name = rdGpArr[1].rstrip().lstrip()
            
    
        if contigName != contig2Name and contig2Start < contig2End  and min(contigStart,contig2Start) < thres and min(lengthDic[contigName]- contigEnd,  lengthDic[contig2Name] - contig2End) < thres:
            conditionForMatch = True
        else:
            conditionForMatch = False

        if conditionForMatch :
            if contigStart < thres:
                dataSet.append((matchLen, contig2Name, contigName))
            
        
        tmp = f.readline()
        
    f.close()  
    dataSet.sort(reverse=True)
    
    print "dataSet:", dataSet[0]
    numberOfContig = len(lengthDic)
    
    leftConnect = [-1 for i in range(numberOfContig)]
    rightConnect = [-1 for i in range(numberOfContig)]
    
    for eachitem in dataSet:
        prefix = eachitem[1].split('_')
        suffix = eachitem[2].split('_')
        
        if prefix[1] == 'p':
            prefixContig = int(prefix[0][6:])*2 
        else:
            prefixContig = int(prefix[0][6:])*2 +1
        
        if suffix[1] == 'p':
            suffixContig = int(suffix[0][6:])*2 
        else:
            suffixContig = int(suffix[0][6:])*2 +1
            
        
        if leftConnect[suffixContig] == -1 and rightConnect[prefixContig] == -1:
            leftConnect[suffixContig] = prefixContig 
            rightConnect[prefixContig] = suffixContig
    
    startList= []
    print leftConnect
    for i in range(len(leftConnect)):
        if leftConnect[i] == -1:
            startList.append(i)
    
    contigList = []
    print startList
    for eachitem in startList:
        tmp = eachitem
        myList = [tmp]
        while rightConnect[tmp] != -1:
            tmp = rightConnect[tmp]
            if tmp != -1:
                myList.append(tmp)
        contigList.append(myList)
            
    i = 0
    fOriginal = open(folderName + "contigs.fasta", 'r')
    readSet = []
    tmp = fOriginal.readline().rstrip()
    tmpRead = ""
    while len(tmp) > 0:
        if tmp[0] == '>':
            if len(tmpRead) >0:
                readSet.append(tmpRead)
                tmpRead = ""
        else:
            tmpRead = tmpRead+ tmp
            
        tmp = fOriginal.readline().rstrip()
    readSet.append(tmpRead)  
    
    seqToPrint = []
    contigUsed = [False for i in range(numberOfContig/2)]
    fAllInOne = open(folderName+"allInOne.fasta",'w' )
    finalList = []
    for eachContig, i in zip(contigList, range(len(contigList))):
        tmpList = []
        for eachitem in eachContig:

            readNum = eachitem/2
            if contigUsed[readNum] ==False:
                seqToPrint.append(eachitem)
                tmpList.append(eachitem)
                contigUsed[readNum] = True
        if len(tmpList) > 0:
            finalList.append(tmpList)
        
    
    fImproved = open(folderName +"improved.fasta", 'w')
    for eachcontig, dummyIndex in zip(finalList, range(len(finalList))):
        fImproved.write(">Segkk"+str(dummyIndex)+'\n')
        for eachseg in eachcontig:
            readNum = eachseg/2
            orientation = eachseg%2

            if orientation == 0:
                fImproved.write(readSet[readNum])
            else:
                fImproved.write(reverseComplement(readSet[readNum]))
        fImproved.write('\n')
        
    fImproved.close()
    
def mainFlow(folderName,mummerLink ):
    greedyAlg(mummerLink, folderName)
    newGraphPipeLine(folderName, mummerLink)
    



t0 = time.time()
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

folderName = sys.argv[1]
mummerLink = sys.argv[2]
mainFlow(folderName,mummerLink)
print  "Time",  time.time() - t0