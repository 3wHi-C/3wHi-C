#!/usr/bin/env python3
from matplotlib import pyplot as plt
import pandas as pd
import os
import json
import sys
if __name__=='__main__':
	filelist=sorted([x for x in os.listdir() if os.path.isdir(x) and x[0]!='.'])
	if len(sys.argv)==1:
		reps=['.']
	else:
		reps=['.']+sys.argv[1:]
	plt.rcParams.update({'font.size':12})
	c={}

	for i in filelist:
		lf=0
		ls=0
		for j in reps:
			with open('%s/%s/fastp.json'%(j,i),'r') as f:
				fpjs=json.load(f)
				lf+=fpjs['summary']['before_filtering']['total_reads']
		with open('%s/%s_hg38XX.bwt2pairs_interaction.sam'%(i,i)) as f:
			ls+=len(f.readlines())
		c[i]=lf/ls
	fig,ax=plt.subplots(figsize=[7,7])
	for i in filelist:
		a=pd.read_table('%s/match_count.txt'%i,header=None)
		plt.plot(a[0][:-1]*c[i],a[2][:-1],label=i)
	#plt.legend(fontsize=10)
	ax.set_xlabel('Rawdata reads number')
	ax.set_ylabel('3-way Complexes number')
	fig.savefig('complex-reads.svg')
	plt.figure().clear()
	fig,ax=plt.subplots(figsize=[7,7])
	for i in filelist:
		a=pd.read_table('%s/match_count.txt'%i,header=None)
		plt.plot(a[0][:-1]*c[i]/3333333,a[2][:-1],label=i)
	#plt.legend(fontsize=10)
	ax.set_xlabel('G Rawdata')
	ax.set_ylabel('3-way Complexes number')
	fig.savefig('complex-Graw.svg')
	plt.figure().clear()
	fig,ax=plt.subplots(figsize=[7,7])
	for i in filelist:
		a=pd.read_table('%s/match_count.txt'%i,header=None)
		plt.plot(a[0][:-1]*c[i],a[2][:-1]/a[0][:-1],label=i)
	ax.set_xlabel('Rawdata reads number')
	ax.set_ylabel('Contribution of rawreads [complexes / reads]')
	#plt.legend(fontsize=10)
	fig.savefig('delta-complex-reads.svg')
	plt.figure().clear()
	fig,ax=plt.subplots(figsize=[7,7])
	for i in filelist:
		a=pd.read_table('%s/match_count.txt'%i,header=None)
		plt.plot(a[0][:-1]*c[i]/3333333,a[2][:-1]/a[0][:-1],label=i)
	ax.set_xlabel('G Rawdata')
	ax.set_ylabel('Contribution of rawreads [complexes / reads]')
	#plt.legend(fontsize=10)
	fig.savefig('delta-complex-Graw.svg')
	plt.figure().clear()
