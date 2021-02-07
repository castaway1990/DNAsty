#!/usr/bin/env python

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
import pysam
import sys
from collections import defaultdict
import numpy as np
import pandas as pd
import itertools
from scipy.sparse import csr_matrix
from pandas.api.types import CategoricalDtype
from scipy.sparse.linalg import svds





def ReadCounter(chunkIdx, bamFile, barcodeList, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
	bam=pysam.AlignmentFile(bamFile, "rb")
	BarcodeSet=set(barcodeList)
	readList=[]
	for indexPos in chunkIdx:
		locus = GenotypeChunkList[indexPos]
		redBases = set()
		readListLocus = []
		for read in bam.fetch(locus[0], locus[1]-1, locus[1]):
		#Check for position coverage
			try:
				CB=read.get_tag("CB")
			except:
				continue
			if not set([CB]).intersection(BarcodeSet):
				continue
			try:
				position = read.positions.index(int(locus[1]) - 1)
			except:
				continue
			try:
				UB=read.get_tag("UB")
			except:
				continue
			if read.query_alignment_qualities[position] < baseQuality:
				continue
			if read.is_secondary:
				continue
			if read.mapq < mapQuality:
				continue
			if read.rlen < readLength:
				continue
			if read.flag == bannedFlag:
				continue
			redBase = read.query_alignment_sequence[position]
			redBases.add(redBase)
			# Max 2 alleles spotted otherwise move to next locus (Biallelic Only)
			if len(redBases) > 2:
				readListLocus = []
				break
			readListLocus.append([locus[0]+ "_"+str(locus[1]), CB, redBase])
		#Min 2 alleles spotted otherwise skip locus pileup (Biallelic Only)
		if len(redBases) == 2:
			readList.extend(readListLocus)
	bam.close()
	print("Chunk number "+ str(chunkIdx) + " completed")
	return readList


def ReadCounterMod(chunkIdx, barcodeList, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
	BarcodeSet=set(barcodeList)
	readList=[]
	for indexPos in chunkIdx:
		locus = GenotypeChunkList[indexPos]
		redBases = set()
		readListLocus = []
		for read in bam.fetch(locus[0], locus[1]-1, locus[1], multiple_iterators = True):
		#Check for position coverage
			try:
				CB=read.get_tag("CB")
			except:
				continue
			if not set([CB]).intersection(BarcodeSet):
				continue
			try:
				position = read.positions.index(int(locus[1]) - 1)
			except:
				continue
			try:
				UB=read.get_tag("UB")
			except:
				continue
			if read.query_alignment_qualities[position] < baseQuality:
				continue
			if read.is_secondary:
				continue
			if read.mapq < mapQuality:
				continue
			if read.rlen < readLength:
				continue
			if read.flag == bannedFlag:
				continue
			redBase = read.query_alignment_sequence[position]
			redBases.add(redBase)
			# Max 2 alleles spotted otherwise move to next locus (Biallelic Only)
			if len(redBases) > 2:
				readListLocus = []
				break
			readListLocus.append([locus[0]+ "_"+str(locus[1]), CB, redBase])
		#Min 2 alleles spotted otherwise skip locus pileup (Biallelic Only)
		if len(redBases) == 2:
			readList.extend(readListLocus)
	bam.close()
	print("Chunk number "+ str(chunkIdx) + " completed")
	return readList


def ReadCounterOverHeadTest(chunkIdx, bamFile, barcodeList, bannedFlag=3844, mapQuality=2, baseQuality=20 , readLength=30):
	BarcodeSet=set(barcodeList)
	readList=[]
	for indexPos in chunkIdx:
		locus = GenotypeChunkList[indexPos]
		time.sleep(2)
	print("Chunk number "+ str(chunkIdx) + " completed")
	return readList


def DFMaker(readList, coverageThreshold, goodBarcodes, CleanGenes):
	Reads=pd.DataFrame(readList, columns = ["Pos","Barcode","Base"])
	#Selecting only locus-barcode minimal coverage entries
	Reads = Reads.groupby(["Pos","Barcode"]).filter(lambda x: len(x) > coverageThreshold)
	print("Defining Alelic combinations per locus...")
	#Defining Al combinations found in dataset
	ReadsStack=Reads[["Pos","Base"]].groupby(["Pos","Base"]).size().unstack(1)
	for cmb in list(itertools.combinations(ReadsStack.columns,2)):
		cmbAlleles = ReadsStack.loc[(ReadsStack[cmb[0]].notna()) & (ReadsStack[cmb[1]].notna())].index
		Reads.loc[Reads.Pos.isin(cmbAlleles),"Al1"] =  cmb[0]
		Reads.loc[Reads.Pos.isin(cmbAlleles),"Al2"] =  cmb[1]
	####LOCI NON UNIQUES CHECK WHY!!!
	####LOCI NON UNIQUES CHECK WHY!!!
	####LOCI NON UNIQUES CHECK WHY!!!
	uniqueLoci = list(set(Reads.Pos))
	####LOCI NON UNIQUES CHECK WHY!!!
	####LOCI NON UNIQUES CHECK WHY!!!
	####LOCI NON UNIQUES CHECK WHY!!!
	###CHECK THIS!! MutReads[["Pos","Base"]].groupby(["Pos","Base"]).size().unstack()
	#Storing Sparse categories
	Pos_c = CategoricalDtype(sorted(uniqueLoci), ordered=True)
	Barcode_c = CategoricalDtype(sorted(list(goodBarcodes)), ordered=True)
	print("Storing Al1Reads...")
	Al1Reads=Reads[Reads["Base"] == Reads["Al1"]][["Pos","Barcode"]]
	Al1Reads=Al1Reads.groupby(["Pos","Barcode"], as_index = False)
	Al1Reads=Al1Reads.apply(len)
	Al1Reads.columns = ["Pos","Barcode","Counts"]
	Al1Reads.index = list(Al1Reads.Pos)
	Pos_Al1 = Al1Reads.Pos.astype(Pos_c).cat.codes
	Barcode_Al1 = Al1Reads.Barcode.astype(Barcode_c).cat.codes
	sparse_Al1 = csr_matrix((Al1Reads["Counts"], (Pos_Al1, Barcode_Al1)), shape=(len(uniqueLoci), len(goodBarcodes)))
	print("Al1Reads stored")
	print("Storing Al2Reads...")
	Al2Reads=Reads[Reads["Base"] == Reads["Al2"]][["Pos","Barcode"]]
	Al2Reads=Al2Reads.groupby(["Pos","Barcode"], as_index = False)
	Al2Reads=Al2Reads.apply(len)
	Al2Reads.columns = ["Pos","Barcode","Counts"]
	Al2Reads.index = list(Al2Reads.Pos)
	Pos_Al2 = Al2Reads.Pos.astype(Pos_c).cat.codes
	Barcode_Al2 = Al2Reads.Barcode.astype(Barcode_c).cat.codes
	sparse_Al2 = csr_matrix((Al2Reads["Counts"], (Pos_Al2, Barcode_Al2)), shape=(len(uniqueLoci), len(goodBarcodes)))
	print("Al1Reads and Al2Reads stored")
	sparseSum=sparse_Al1+sparse_Al2
	sparseSum.data = 1/sparseSum.data
	SparseBval=sparse_Al1.multiply(sparseSum)
	#SparseBval[np.isnan(SparseBval)] = 0
	#SparseBval=csr_matrix(SparseBval)
	return SparseBval
