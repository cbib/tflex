import os


# @author johaGL

# 'strucfile' option must be specified in configfile:
#     note :  a {pair name} example is 'R2' or '2' or equivalent
#   1 -- > if {sampleid}_{pair name}.fastq.gz
#         examples: D123546_1.fastq.gz , T9_876-sansglu_R1.fastq.gz
#   2 -- > if {sampleid}_{pair name}_{run}.fastq.gz, and {run}
#         example: T9_876-sansglu_R1_001.fastq.gz , so {run} is '001'

def makelistsunique(FASTQNAMES, PAIRNAMES, RL):
    # FASTQNAMES lists all file names without pair name (== twice sample name). Make it unique SAMPLES:
    SAMPLES = list(set(FASTQNAMES)) # ok unique so we have true sample names in SAMPLES list
    # PAIRNAMES lists n/2 times each pairname, make unique and order into PAIRS
    PAIRS = sorted(set(PAIRNAMES)) # expected result is [ 1,  2 ] or ['R1' , 'R2'] or equivalent
    RUN = sorted(set(RL))
    print(f"  --> detected a number of {len(SAMPLES)} samples (paired-end)")
    print(f"  --> printing a sample name :  {SAMPLES[0]} ")
    print(f"  --> pair names are :  {PAIRS}")
    return SAMPLES, PAIRS, RUN

def listOneelemTostrdesired(strucfile, RUN):
    """ requires list to be lenght 1. Yields str preceded by '_' if strucfile is 2"""
    assert len(set(RUN)) == 1, "more than one element, this function NEEDS list length 1"
    if strucfile == 1:
        return RUN[0] # does nothing, kept empty
    elif strucfile == 2:
         return "_"+RUN[0]
