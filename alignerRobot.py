import os
import sys
import houseKeeper
from multiprocessing import Pool
import time
#For SLURM task array blocking until completion (detected by sentinel file):
import fcntl

def extractMumData(folderName, fileName):
    # "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    f = open(folderName + fileName, 'r')
    dataList = []
    
    for i in range(6):
        tmp = f.readline()

    while len(tmp) > 0:
        info = tmp.split('|')
        filterArr = info[1].split()
        rdGpArr = info[-1].split('\t')
        firstArr = info[0].split()
        
        matchLenArr = info[2].split()
    
        matchLen1 = int(matchLenArr[0])
        matchLen2 = int(matchLenArr[1])    
        percentMatch = float(info[3])
        
        
        helperStart, helperEnd = int(firstArr[0]), int(firstArr[1])
        readStart, readEnd = int(filterArr[0]) , int(filterArr[1])
        
        helperName = rdGpArr[0].rstrip().lstrip()
        readName = rdGpArr[1].rstrip().lstrip()
        
        dataList.append([helperStart, helperEnd , readStart, readEnd, matchLen1, matchLen2, percentMatch, helperName, readName ])
    
        tmp = f.readline().rstrip()
                
    f.close()
    
    return dataList



def useMummerAlign(mummerLink, folderName, outputName, referenceName, queryName, specialForRaw = False, specialName = "", refinedVersion= False):
    nucmerMummer(specialForRaw, mummerLink, folderName, outputName, referenceName, queryName, refinedVersion)
    showCoorMummer(specialForRaw, mummerLink, folderName, outputName, specialName)
    
    

def nucmerMummer(specialForRaw, mummerLink, folderName, outputName, referenceName, queryName,refinedVersion):
    if not refinedVersion:
        if not specialForRaw:
            if houseKeeper.globalFast:
                command = mummerLink + "nucmer -b 50  --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " + folderName + queryName
            else:
                command = mummerLink + "nucmer --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " + folderName + queryName
        else:
            if houseKeeper.globalFast:
                command = mummerLink + "nucmer -b 50  --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " +  queryName
            else:
                command = mummerLink + "nucmer --maxmatch -p " + folderName + outputName + " " + folderName + referenceName + " " +  queryName
    else:
        command = mummerLink + "nucmer -l 10 --maxmatch -p " + folderName + outputName +" " + folderName +referenceName +" " + queryName
        
    os.system(command)
    

def showCoorMummer(specialForRaw, mummerLink, folderName, outputName, specialName):
    if not specialForRaw:
        command = mummerLink + "show-coords -r " + folderName + outputName + ".delta > " + folderName + outputName + "Out"
    else:
        command = mummerLink + "show-coords -r " + folderName + outputName + ".delta > " + folderName + specialName 
        
    os.system(command)


#Defunct method?
def combineMultipleCoorMum(specialForRaw, mummerLink, folderName, outputName, specialName, numberOfFiles):
    print ""
    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
        showCoorMummer(specialForRaw, mummerLink, folderName, outputName+indexOfMum, specialName+indexOfMum)
    
    
    command =  "head -5 "+ folderName + specialName +"01" + "> " + folderName +specialName
    os.system(command)

    for dummyI in range(1, numberOfFiles + 1):
        indexOfMum = ""
        if dummyI < 10:
            indexOfMum = "0" + str(dummyI)
        else:
            indexOfMum = str(dummyI)
    
        command = " tail -n+6 "+ folderName + specialName +indexOfMum +">> " + folderName +specialName
        os.system(command)

def zeropadding(i):
    tmpi = ""

    if i < 10:
        tmpi = "0" + str(i)
    else:
        tmpi = str(i)
    return tmpi
   
   
def useMummerAlignBatch(mummerLink, folderName, workerList, nProc ,specialForRaw = False, refinedVersion = False):
    # Format for workerList : [[outputName, referenceName, queryName, specialName]... ]
    # nProc : a parameter on how many threads should be created each time
    # Goal : parallelize this part  
    
    
    if not houseKeeper.globalLarge:
        p = Pool(processes=nProc)
        results = []
        
        for eachitem in workerList:
            outputName, referenceName, queryName, specialName = eachitem
            results.append(p.apply_async(useMummerAlign, args=(mummerLink, folderName, outputName, referenceName, queryName, specialForRaw , specialName, refinedVersion)))
        
        outputlist = [itemkk.get() for itemkk in results]
        print  len(outputlist)
        p.close()
    else:
        '''
        a) Split
        b) align several times
        c) join the query
        raw_reads.part-01
        '''
        p = Pool(nProc)
        results = []
        numberRefFiles = houseKeeper.globalNumRefFiles
        numberQueryFiles = houseKeeper.globalNumQueryFiles
#        numberOfFiles = 10
        
        for eachitem in workerList:
            print eachitem
            outputName, referenceName, queryName, specialName = eachitem[0], eachitem[1], eachitem[2] , eachitem[3]
            #Workaround for MUMmering a file against itself when splitting reference and query into different
            #numbers of parts:
            if referenceName == queryName: #If MUMming a file against itself
                #Make a symlink to the input file to act as a query file for splitting:
                symQueryName = queryName.replace(".fasta", "-symlink.fasta")
                os.symlink(queryName, symQueryName)
                queryName = symQueryName
            
            bindir =  os.path.abspath(os.path.dirname(sys.argv[0]))   
            command = bindir + "/fasta-splitter.pl --n-parts " + str(numberRefFiles) + " " + folderName + referenceName
            os.system(command)
            

            if specialForRaw : 
                queryNameMod = queryName
            else:
                queryNameMod = folderName + queryName

            command = bindir + "/fasta-splitter.pl --n-parts " + str(numberQueryFiles) + " " + queryNameMod
            os.system(command)
            
        if houseKeeper.globalUseSlurm == True:
            os.mkdir('locks')
            locks_parent_directory = os.path.abspath(os.getcwd())
            print locks_parent_directory+'/locks/'
        
            
        for eachitem in workerList:   
            outputName, referenceName, queryName, specialName = eachitem[0], eachitem[1], eachitem[2] , eachitem[3]
            #Workaround for MUMmering a file against itself when splitting reference and query into different
            #numbers of parts:
            if referenceName == queryName: #If MUMming a file against itself
                #Use the symlink to the input file:
                queryName = queryName.replace(".fasta", "-symlink.fasta")
                
            if houseKeeper.globalUseSlurm == True:
                if houseKeeper.globalFast:
                    nucmer_fast_option = "-b 50 "
                else:
                    nucmer_fast_option = ""
                    
                if refinedVersion:
                    nucmer_fast_option = ""
                    nucmer_refined_option = "-l 10 "
                else:
                    nucmer_refined_option = ""
                    
                slurmscript = open(outputName+'_mummer_slurm.sh', 'w')
                slurmscript.write('#!/bin/bash\n')
                slurmscript.write('JOBID=$SLURM_ARRAY_JOB_ID\n')
                slurmscript.write('ID=$SLURM_ARRAY_TASK_ID\n')
                slurmscript.write('PADDEDID=$(printf %02d ${ID})\n')
                slurmscript.write('pwd\n')
                slurmscript.write('exec 7>'+locks_parent_directory+'/locks/${JOBID}_${ID}\n')
                slurmscript.write('flock 7\n')
                slurmscript.write('REF="'+referenceName[0:-6]+'.part-${PADDEDID}.fasta"\n')
                slurmscript.write('for QUERYID in {01..'+str(numberQueryFiles)+'}\n')
                slurmscript.write('   do ')
                if specialForRaw:
                    slurmscript.write('QUERY="'+queryName[0:-6]+'-${QUERYID}.fasta"\n')
                else:
                    slurmscript.write('QUERY="'+queryName[0:-6]+'.part-${QUERYID}.fasta"\n')
                    
                slurmscript.write('   '+mummerLink+'nucmer '+nucmer_refined_option+nucmer_fast_option+'--maxmatch -p '+folderName+outputName+'${PADDEDID}${QUERYID} ${REF} ${QUERY}')
                slurmscript.write('done\n')
                slurmscript.write('flock -u 7\n')
                slurmscript.close()
                command = 'sbatch --array=1-'+str(numberRefFiles)+houseKeeper.globalSlurmParams+' '+outputName+'_mummer_slurm.sh'
                os.system(command)
            else:
                for i in range(1, numberRefFiles+1):
                    for j in range(1, numberQueryFiles+1):
                        if specialForRaw : 
                            tmpRefName , tmpQryName = referenceName[0:-6] + ".part-" + zeropadding(i) +".fasta",  queryName[0:-6] + "-" + zeropadding(j) + ".fasta"
                        else:
                            tmpRefName , tmpQryName = referenceName[0:-6] + ".part-" + zeropadding(i) +".fasta",  queryName[0:-6] + ".part-" + zeropadding(j) + ".fasta"
                    
                        results.append(p.apply_async(nucmerMummer, args =(specialForRaw, mummerLink, "", folderName + outputName +zeropadding(i)+zeropadding(j), tmpRefName, tmpQryName, refinedVersion)))
      
      
        
        if houseKeeper.globalUseSlurm == True:
            #Wait for the sentinel files to release their locks:
            num_done = 0
            lockfiles = os.listdir(locks_parent_directory+'/locks')
            for lockfile in lockfiles:
                current_lock = open(locks_parent_directory+'/locks/'+lockfile, 'rb')
                flock(current_lock, fcntl.LOCK_EX)
                flock(current_lock, fcntl.LOCK_UN)
                current_lock.close()
                os.unlink(locks_parent_directory+'/locks/'+lockfile)
                num_done += 1
            
            os.rmdir(locks_parent_directory+'/locks')
            print num_done
        else:
            outputlist = [itemkk.get() for itemkk in results]
            print len(outputlist)
        
        p.close()
        
        for eachitem in workerList:
            outputName, referenceName, queryName, specialName = eachitem           
            if not specialForRaw:
                outNameMod =  folderName + outputName + "Out"
            else:
                outNameMod = folderName + specialName 
        
        
            tmpName = folderName + outputName +zeropadding(1)+zeropadding(1) + ".delta"
        
            command = mummerLink + "show-coords -r " + tmpName + "| head -5 > " + outNameMod
            os.system(command)
            
            for i in range(1, numberRefFiles+1):
                for j in range(1, numberQueryFiles+1):
                    
                    tmpName = folderName + outputName +zeropadding(i)+zeropadding(j) + ".delta"
                    command = mummerLink + "show-coords -r " + tmpName + "| tail -n+6 >> " + outNameMod
                    os.system(command)
      #Should we be deleting all the split files after collating all the delta files together?
      #It might be a good idea to reduce the clutter in the working directory.
      #Perhaps also delete the partial delta files after successful concatenation.
        

def transformCoor(dataList):
    # "Format of the dataList :  1      765  |    11596    10822  |      765      775  |    84.25  | ref_NC_001133_       scf7180000000702"
    newList = []
    
    for eachitem in dataList:
        if eachitem[2] < eachitem[3]:
            newList.append(eachitem)
        else:
            tmpitem = eachitem
            tmp = tmpitem[2]
            tmpitem[2] = tmpitem[3]
            tmpitem[3] = tmp
            newList.append(tmpitem)
    
    return newList
