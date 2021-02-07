#!/usr/bin/env python

import argparse
import os
import scanpy as sc
import pysam
import pandas as pd
import time
import itertools
import io
import numpy as np
from multiprocessing import Pool
import sys
from collections import defaultdict
import re
import math
from collections import Counter
from pybedtools import BedTool






bamFile = "/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20273_158/outs/possorted_genome_bam.bam"
filtered_feature_bc_matrix = "/hpcnfs/scratch/temporary/Dav_vc/0_CellrangerAlign/Sample_S20273_158/outs/filtered_feature_bc_matrix/"
IDsPath = "/hpcnfs/scratch/temporary/Dav_vc/20_Final_DMX/Sample_S20273_158/SCanSNP/doubletsMarked.tsv"
annotationGTF = "/hpcnfs/techunits/bioinformatics/refdata/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf"
commonVarsPath = "/hpcnfs/scratch/GT/References/Ensembl/1000G/ALL.wgs.shapeit2_integrated_snvindels_v2a.GRCh38_COMMON.bed2"
nThreads = 10


HK_genes,barcodeList = ExtractHKgenes(filtered_feature_bc_matrix, IDsPath, SelectedID="MIFF1")
CleanGenes = ExtractLociList(HK_genes, annotationGTF, commonVarsPath)


GenotypeChunkDict=ChunkMaker(CleanGenes, nThreads)
#Converting genotypeschunk into list
GenotypeChunkList,GenotypeChunkIndexesList = FlattenDict(GenotypeChunkDict)

del GenotypeChunkDict

#FireUp main Pileupper

start_time = time.time()

results = []
pool=Pool(nThreads)

for chunkIdx in GenotypeChunkIndexesList:
	result = pool.apply_async(ReadCounter, (chunkIdx, bamFile, barcodeList))
	results.append(result)


pool.close()
pool.join()

print("--- %s seconds ---" % (time.time() - start_time))



Reads=list(itertools.chain.from_iterable(([result.get() for result in results])))

# with open('/hpcnfs/scratch/temporary/ieo4777/DNAsty/ReadsFile.p', 'wb') as fp:
#     pickle.dump(Reads, fp)
#
# with open ('/hpcnfs/scratch/temporary/ieo4777/DNAsty/ReadsFile.p', 'rb') as fp:
#     test = pickle.load(fp)



test=DFMaker(Reads, 10, barcodeList, CleanGenes)
