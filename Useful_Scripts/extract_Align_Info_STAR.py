#!usr/bin/env python

##
## Parse the aligned samples-folders and extract info for each bam:
## name, size, mapping from _Log.final.out
##

import os
import sys
import numpy as np
import pandas as pd
import fnmatch

workdir="/data/Maria/Mbnl1KO_Xaviere/Aligned/"
#os.chdir(workdir)
sampledirs = [s for s in os.listdir(workdir)]
file_info = pd.DataFrame(columns=['bam_Size', 'Uniq_Mapped_pct'],index=sampledirs)
## Read bam file for each sample:
for s in sampledirs :
    "\nSample {s} being read...\n"
    f_in_dir = os.listdir(workdir+s)
    suffix="*_Aligned.sortedByCoord.out.bam"
    #extract the bam from subdir
    found = [fnmatch.fnmatch(f, suffix) for f in os.listdir(workdir+s)]
    i = found.index(True)
    f"\nBam file found in {i}-th position.\n"
    #a/ show size in Gb
    mybam = os.listdir(workdir+s)[i]
    fstats = os.stat(workdir+s+"/"+mybam)
    file_size = round(fstats.st_size / 1024 **3, 2 )
    file_info.loc[s]["bam_Size"] = file_size
    #b/ Read from "_Log.final.out" --> "Uniquely mapped reads %"
    log = "*_Log.final.out"
    found = [fnmatch.fnmatch(f,log) for f in os.listdir(workdir+s)]
    j = found.index(True)
    f"\nLog file found in {j}-th position.\n"
    log_file = os.listdir(workdir+s)[j]
    with open(workdir+s+"/"+log_file) as f:
        for line in f :
            if line.find("Uniquely mapped reads") != -1 : # = match found !
                pct_mapped = line.split("\t")[1].split("\n")[0]
                file_info.loc[s]['Uniq_Mapped_pct'] = pct_mapped

#end.for

file_info.head(5)

## Save to a file :
file_info.to_csv(workdir+"bam_file_info.tsv",sep="\t",encoding='utf-8')
