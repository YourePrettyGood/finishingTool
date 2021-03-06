
import time
import argparse
import os

import nonRedundantResolver
import overlapResolver
import gapFiller
import twoRepeatOneBridgeSolver
import houseKeeper

###################################################### Starting point
def mainFlow(folderName , mummerLink, pickupname, mapcontigsname):      
    print "Go Bears! ! !" 
    
    print "pickupname, mapcontigsname", pickupname, mapcontigsname
    
    
    if not pickupname in ["noEmbed.fasta", "improved.fasta", "improved2.fasta"]:
        print "NonRedundantResolver RemoveEmbedded"
        nonRedundantResolver.removeEmbedded(folderName , mummerLink)
     
    if not pickupname in ["improved.fasta", "improved2.fasta"]:
        print "OverlapResolver FetchSuccessor"
        overlapResolver.fetchSuccessor(folderName , mummerLink)
        print "OverlapResolver FormSeqGraph"
        overlapResolver.formSeqGraph(folderName , mummerLink)
    
    if not pickupname in ["improved2.fasta"]:
        print "GapFiller FillGap"
        gapFiller.fillGap(folderName , mummerLink)
    
    print "TwoRepeatOneBridgeSolver XPhased"
    twoRepeatOneBridgeSolver.xPhased(folderName , mummerLink)
    
    # ECReduction(folderName , mummerLink )
    # compareWithReference(folderName , mummerLink)
    
    if mapcontigsname != None:
        print "HouseKeeper PerformMapping"
        houseKeeper.performMapping(folderName, mummerLink, mapcontigsname)
        
    print "<3 Do cool things that matter <3"

# folderName = "S_cerivisea/"
# mummerLink = "MUMmer3.23/"

t0 = time.time()

parser = argparse.ArgumentParser(description='FinisherSC : a repeat-aware tool to upgrade de-novo assembly with long reads')
parser.add_argument('folderName')
parser.add_argument('mummerLink')
parser.add_argument('-p', '--pickup', help='Picks up existing work (input is noEmbed.fasta, improved.fasta or improved2.fasta)', required=False)
parser.add_argument('-o', '--mapcontigs', help='Maps new contigs to old contigs(input is of the format of contigs.fasta_improved3.fasta which means improved3.fasta will be mapped back to contigs.fasta; Output can be found in mappingResults.txt in the destinedFolder;)', required=False)
parser.add_argument('-f', '--fast', help= 'Fast aligns contigs (input is True)', required=False)
parser.add_argument('-par', '--parallel', help= 'Fast aligns contigs (input is maximum number of threads)', required=False)
parser.add_argument('-l', '--large', help= 'Large number of contigs/large size of contigs (input is True)', required=False)
parser.add_argument('-r', '--numreffiles', help='Split MUMmer reference file into n parts (default is 10)', required=False)
parser.add_argument('-q', '--numqueryfiles', help='Split MUMmer query file into n parts (default is 10)', required=False)
parser.add_argument('-s', '--use_slurm', help='Submit MUMmer jobs as a SLURM task array (default is False)', required=False)
parser.add_argument('-t', '--slurm_timeout', help='Timeout to use for SLURM (default is 24 hours)', required=False)
parser.add_argument('-m', '--slurm_mem', help='Memory to allocate for each task in SLURM array (default is 32000 MB)', required=False)
parser.add_argument('-c', '--slurm_cpus', help='CPUs per task to allocate for each task in SLURM array (default is 1)', required=False)
parser.add_argument('-z', '--slurm_misc', help='Miscellaneous SLURM parameters (no default)', required=False)


args = vars(parser.parse_args())

print "args", args
pathExists, newFolderName, newMummerLink = houseKeeper.checkingPath(args['folderName'], args['mummerLink'])

if args['fast'] == "True":
    houseKeeper.globalFast = True
else:
    houseKeeper.globalFast = False

if args['parallel'] != None:
    houseKeeper.globalParallel = int(args['parallel'])
else:
    houseKeeper.globalParallel = 1


if args['large'] == "True":
    houseKeeper.globalLarge = True
else:
    houseKeeper.globalLarge = False

#Allow for the reference files for MUMmer to be split into a custom number of parts to control memory usage:
#Note: Per Adam Phillippy, the suffix tree built for the reference in MUMmer requires, as a rule of thumb,
# 17 bytes per base
if args['numreffiles'] != None:
    houseKeeper.globalNumRefFiles = int(args['numreffiles'])
else:
    houseKeeper.globalNumRefFiles = 10

if args['numqueryfiles'] != None:
    houseKeeper.globalNumQueryFiles = int(args['numqueryfiles'])
else:
    houseKeeper.globalNumQueryFiles = 10

#To deal with jobs that require tons of MUMmer calls, allow for SLURM task array submission:
if args['use_slurm'] == "True":
    houseKeeper.globalUseSlurm = True
    os.mkdir('locks')
else:
    houseKeeper.globalUseSlurm = False
    
if args['slurm_timeout'] != None:
    houseKeeper.globalSlurmParams += ' -t '+args['slurm_timeout']
else:
    houseKeeper.globalSlurmParams += ' -t 24:00:00'
    
if args['slurm_mem'] != None:
    houseKeeper.globalSlurmParams += ' --mem='+args['slurm_mem']
else:
    houseKeeper.globalSlurmParams += ' --mem=32000'
    
if args['slurm_cpus'] != None:
    houseKeeper.globalSlurmParams += ' --cpus-per-task='+args['slurm_cpus']
else:
    houseKeeper.globalSlurmParams += ' --cpus-per-task=1'
    
if args['slurm_misc'] != None:
    houseKeeper.globalSlurmParams += ' '+args['slurm_misc']


if pathExists:
    print "Input parameters:"
    print "Fast?: " + str(houseKeeper.globalFast)
    print "Number of CPU cores: " + str(houseKeeper.globalParallel)
    print "Large genome?: " + str(houseKeeper.globalLarge)
    print "Split MUMmer reference into how many files?: " + str(houseKeeper.globalNumRefFiles)
    print "Split MUMmer query into how many files?: " + str(houseKeeper.globalNumQueryFiles)
    print "MUMmer path: " + newMummerLink
    print "Working directory: " + newFolderName
    print "Pick up from which step?: " + str(args['pickup'])
    print "Map new contigs to old contigs?: " + str(args['mapcontigs'])
    print "Split raw reads into how many files?: " + str(max(20, houseKeeper.globalParallel))
    print "Use SLURM task array for MUMmer?: " + str(houseKeeper.globalUseSlurm)
    print "SLURM task submission parameters: " + str(houseKeeper.globalSlurmParams)
    mainFlow(newFolderName, newMummerLink, args['pickup'], args['mapcontigs'])
else:
    print "Sorry. The above folders or files are missing. If you continue to have problems, please contact me(Ka-Kit Lam) at kklam@eecs.berkeley.edu"

print  "Time", time.time() - t0
