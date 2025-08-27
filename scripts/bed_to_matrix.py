#!/usr/bin/env python3
import pandas as pd
from collections import Counter
import sys,os
prefix=sys.argv[1]
resolutions=set(map(lambda x:int(x[:x.index('_')]),os.listdir(sys.argv[2])))
a1=pd.read_csv(f'{prefix}_0.bed',sep='\t',header=None)
a2=pd.read_csv(f'{prefix}_1.bed',sep='\t',header=None)
def bed3(resolution):
	bed=pd.read_csv(f'{sys.argv[2]}/{resolution}_abs.bed',sep='\t',header=None)
	b={x:{} for x in bed[0]}
	for _,i in bed.iterrows():
		b[i[0]][i[1]//resolution]=i[3]
	del bed
	c=pd.DataFrame()
	c[0]=a1.apply(lambda x:b[x[0]][x[1]//resolution],axis=1)
	c[1]=a2.apply(lambda x:b[x[0]][x[1]//resolution],axis=1)
	c[2]=a3.apply(lambda x:b[x[0]][x[1]//resolution],axis=1)
	del b
	d=['\t'.join(map(str,sorted(x))) for _,x in c.iterrows()]
	d_2=['\t'.join(map(str,sorted(y))) for _,x in c.iterrows() for y in [[x[0],x[1]],[x[1],x[2]],[x[2],x[0]]]]
	del c
	e=Counter(d)
	del d
	f=sorted([list(map(int,x.split()))+[y] for x,y in e.items()])
	with open(f'{prefix}_w3_{resolution}.matrix3','w') as o:
		o.write(''.join(['\t'.join(map(str,x))+'\n' for x in f]))
	e=Counter(d_2)
	f=sorted([list(map(int,x.split()))+[y] for x,y in e.items()])
	with open(f'{prefix}_w3_{resolution}.matrix','w') as o:
		o.write(''.join(['\t'.join(map(str,x))+'\n' for x in f]))
		
def bed2(resolution):
	bed=pd.read_csv(f'{sys.argv[2]}/{resolution}_abs.bed',sep='\t',header=None)
	b={x:{} for x in bed[0]}
	for _,i in bed.iterrows():
		b[i[0]][i[1]//resolution]=i[3]
	del bed
	c=pd.DataFrame()
	c[0]=a1.apply(lambda x:b[x[0]][x[1]//resolution],axis=1)
	c[1]=a2.apply(lambda x:b[x[0]][x[1]//resolution],axis=1)
	del b
	d=['\t'.join(map(str,sorted(x))) for _,x in c.iterrows()]
	del c
	e=Counter(d)
	del d
	f=sorted([list(map(int,x.split()))+[y] for x,y in e.items()])
	with open(f'{prefix}_w3_{resolution}.matrix','w') as o:
		o.write(''.join(['\t'.join(map(str,x))+'\n' for x in f]))
try:
	a3=pd.read_csv(f'{prefix}_2.bed',sep='\t',header=None)
except:
	bed_to_matrix=bed2
else:
	bed_to_matrix=bed3


if __name__=='__main__':
	from multiprocessing import Pool
	#with Pool(len(resolutions)) as pool:
	#	pool.map(bed_to_matrix,resolutions)
	for i in [5000,10000,25000,10000]:
		bed_to_matrix(i)
