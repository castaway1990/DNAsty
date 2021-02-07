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
import networkx as nx
import scanpy as sc
from pybedtools import BedTool



def ExtractHKgenes(filtered_feature_bc_matrix, IDsPath, SelectedID, MinCellsPCT=.4, TopMeanVar=150):
	'''
	Experiment housekeeping detection module
	'''
	adata = sc.read_10x_mtx(filtered_feature_bc_matrix,var_names='gene_symbols',cache=True)
	IDsDF = pd.read_csv(IDsPath, sep="\t", index_col=0)
	adata.obs['ID'] = IDsDF.loc[adata.obs.index, "ID_Qual"]
	adata = adata[adata.obs.ID == SelectedID]
	# Min Cells
	sc.pp.filter_cells(adata, min_counts=2000)
	# Min Cells
	minCellsHouseKeeping=round(len(adata.obs)*MinCellsPCT)
	sc.pp.filter_genes(adata, min_cells=minCellsHouseKeeping)
	#Top N genes mean/variance
	mean=sc.pp._utils._get_mean_var(adata.X)[0]
	varIance=sc.pp._utils._get_mean_var(adata.X)[1]
	TopPos=np.argsort(mean/varIance)[-(TopMeanVar):]
	adata=adata[:,TopPos]
	#Return OBS and VARs
	HK_genes=pd.Series(adata.var.index)
	goodBarcodes = adata.obs_names
	return HK_genes,goodBarcodes





def ExtractLociList(HK_genes, annotationGTF, commonVarsPath):
	'''
	we extract loci list starting from HKgenes, gtf purged from dbSNP common variants
	'''
	with open(annotationGTF, 'r') as f:
		header = 0
		for line in f:
			if line.startswith('#'):
				header = header+1
			else:
				break
	genes=pd.read_csv(annotationGTF, sep = "\t", skiprows = header, header=None, usecols = [0,2,3,4,8], names = ["CHROM","FEATURE","POS","END","ANNO"])
	genes=genes[genes["FEATURE"] == "transcript"]
	genes["ANNO"]=genes["ANNO"].str.split('gene_name').str[1].str.split(";").str[0].str.replace('"', '').str.replace(' ', '')
	genes = genes[genes["ANNO"].isin(list(HK_genes))]
	del genes["FEATURE"]
	genesBed=BedTool.from_dataframe(genes)
	genesBed=genesBed.sort().merge()
	commonVarsBed = BedTool(commonVarsPath)
	CleanGenes=genesBed.subtract(commonVarsBed).to_dataframe()
	CleanGenes.columns = ["CHROM","POS","END"]
	CleanGenes["len"] = -(CleanGenes.POS - CleanGenes.END)
	#CleanGenes["FragID"] = CleanGenes.start.astype(str)+"_"+CleanGenes.end.astype(str)
	Range = [item for sublist in [list(range(length)) for length in CleanGenes.len] for item in sublist]
	CleanGenes=CleanGenes.loc[CleanGenes.index.repeat(CleanGenes.len)]
	CleanGenes.len = Range
	del CleanGenes["END"]
	CleanGenes.POS = CleanGenes.POS - CleanGenes.len
	del CleanGenes["len"]
	CleanGenes.index = CleanGenes.CHROM.astype(str)+"_"+CleanGenes.POS.astype(str)
	####LOCI NON UNIQUES CHECK WHY!!!
	####LOCI NON UNIQUES CHECK WHY!!!
	####LOCI NON UNIQUES CHECK WHY!!!
	CleanGenes = CleanGenes[~CleanGenes.index.duplicated(keep='first')]
	####LOCI NON UNIQUES CHECK WHY!!!
	####LOCI NON UNIQUES CHECK WHY!!!
	####LOCI NON UNIQUES CHECK WHY!!!
	return CleanGenes




def splitContig(contig, genotypes, nThreads):
	ChunkDICT = {}
	RemainsDF = pd.DataFrame()
	GenotypeChunk = genotypes[genotypes["CHROM"] == contig].sample(frac = 1)
	if len(GenotypeChunk) >= nThreads:
		chunkSizes = math.floor(len(GenotypeChunk)/int(nThreads))
		LeftOver = len(GenotypeChunk)%nThreads
		ChunkStart = 0
		for Chunk in range(0,nThreads):
			if Chunk < nThreads-1:
				ChunkDICT[Chunk] = GenotypeChunk.iloc[range(ChunkStart, ChunkStart+chunkSizes)]
			elif Chunk == nThreads-1:
				ChunkDICT[Chunk] = GenotypeChunk.iloc[range(ChunkStart, ChunkStart+chunkSizes+LeftOver)]
			ChunkStart = ChunkStart + chunkSizes
	else:
		RemainsDF = RemainsDF.append(GenotypeChunk)
	return ChunkDICT, RemainsDF





def ChunkMaker(CleanGenes, nThreads):
	'''
	Divide loci list into chunk for parallel processing
	'''
	GenotypeChunkDICT = {}
	Remains = pd.DataFrame()
	GenotypeChunkDICT_final = {}
	CleanGenes_sampled = CleanGenes.sample(frac = 1)
	for contig in list(CleanGenes["CHROM"].unique()):
		GenotypeChunkDICT[contig] = splitContig(contig, CleanGenes_sampled, nThreads)[0]
		Remains = Remains.append(splitContig(contig, CleanGenes_sampled, nThreads)[1])
	for Chunk in range(0,nThreads):
		GenotypeChunkDICT_final[Chunk] = pd.DataFrame()
		for key in GenotypeChunkDICT.keys():
			try:
				GenotypeChunkDICT_final[Chunk] = GenotypeChunkDICT_final[Chunk].append(GenotypeChunkDICT[key][Chunk])
			except:
				continue
	minChunk=min(GenotypeChunkDICT_final, key=lambda k: len(GenotypeChunkDICT_final[k]))
	GenotypeChunkDICT_final[minChunk] = GenotypeChunkDICT_final[minChunk].append(Remains)
	return GenotypeChunkDICT_final



def FlattenDict(GenotypeChunkDict):
	'''
	Transform dictionary into List of lists + ranges for efficiency
	'''
	GenotypeChunkList = []
	GenotypeChunkIndexesList = []
	ChunkStart = 0
	for k in GenotypeChunkDict.keys():
		chunk=GenotypeChunkDict[k].loc[:,["CHROM","POS"]]
		# Create Chunk Index and append to list
		GenotypeChunkIndexesList.append(range(ChunkStart, ChunkStart+chunk.shape[0]))
		ChunkStart=ChunkStart+chunk.shape[0]
		# Re-format chunk to list-like object
		chunk["POS"] = chunk["POS"].astype(int)
		chunk=chunk.values.tolist()
		GenotypeChunkList.append(chunk)
	GenotypeChunkList = [item for sublist in GenotypeChunkList for item in sublist]
	return GenotypeChunkList, GenotypeChunkIndexesList
