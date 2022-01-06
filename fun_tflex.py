import os


# @author johaGL

# 'strucfile' option must be specified in configfile:
#     note :  a {pair name} example is 'R2' or '2' or equivalent
#   1 -- > if {sampleid}_{pair name}.fastq.gz
#         examples: D123546_1.fastq.gz , T9_876-sansglu_R1.fastq.gz
#   2 -- > if {sampleid}_{pair name}_{run}.fastq.gz, and {run}
#         example: T9_876-sansglu_R1_001.fastq.gz , so {run} is '001'

def makelistsunique(FASTQNAMES, PAIRNAMES, RL): 
    """prints check-in messages, outputs 3 lists: SAMPLES, PAIRS, RUN  """
    # FASTQNAMES lists all file names without pair name (== twice sample name). Make it unique SAMPLES:
    SAMPLES = list(set(FASTQNAMES)) # ok unique so we have true sample names in SAMPLES list
    # PAIRNAMES lists n/2 times each pairname, make unique and order into PAIRS
    PAIRS = sorted(set(PAIRNAMES)) # expected result is [ 1,  2 ] or ['R1' , 'R2'] or equivalent
    RUN = sorted(set(RL))
    #print(f"  --> detected a number of {len(SAMPLES)} samples (paired-end)")
    #print(f"  --> printing a sample name :  {SAMPLES[0]} ")
    #print(f"  --> pair names are :  {PAIRS}")
    return SAMPLES, PAIRS, RUN

def listtostr_single(strucfile, RUN):
    """ requires list to be length 1. Yields str preceded by '_' if strucfile is 2"""
    assert len(set(RUN)) == 1, "more than one element, this function NEEDS list length 1"
    if strucfile == 1:
        return RUN[0] # does nothing, kept empty
    elif strucfile == 2:
         return "_"+RUN[0] 

def d2o(afilename, alistspecies):
    """
    input : csv file. Delimiter MUST BE comma
    output : Dictionary, key is species, pointing to samples list      
    """
    delim=","
    ds = {}
    with open(afilename, 'r') as f:
        texto = f.readlines()
    for l in texto:
        lol = l.strip().split(delim)
        samplename = lol[0]
        speciesname = lol[1]
        if speciesname not in ds.keys():
            ds[speciesname] = [samplename]
        else:
            ds[speciesname].append(samplename)
    assert set(ds.keys()) == set(alistspecies), \
          "ERROR, csv file contains species other than those in \
            config file 'species' key (maybe spelling error?)"
    return ds

def dicosamples_species(afilename, alistspecies):
    """input: csv file. Delimiter MUST BE comma
       output : Dictionary, samples are keys pointing to species they belong to
    """
    delim = ","
    ds = {}
    with open(afilename,'r') as f:
        texto = f.readlines()
    detectedspecies = []   
    for l in texto:
        lol = l.strip().split(delim)
        samplename = lol[0]
        speciesname = lol[1]
        if samplename not in ds.keys():
            ds[samplename] = {'sp' : speciesname}
            detectedspecies.append(speciesname)
    assert set(detectedspecies) == set(alistspecies), \
          "ERROR, csv file contains species other than those in \
           config file 'species' field (spelling error?)"
    return ds

