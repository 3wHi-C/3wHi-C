#!/usr/bin/env python3
import numpy as np
import pandas as pd
import sys
chros=['chr%s'%x for x in range(1,23)]+['chrX']
if __name__=="__main__":
	matrix3=pd.read_table(sys.argv[1],header=None)
	b=pd.read_table(sys.argv[2],header=None)
	b=b[b[0].apply(lambda x:x in chrm)][3]
	binmax=b[3].max()
	mp=dict(zip(b[3],b[0]))
	MP={x:[] for x in mp.values()}
	matrix3=matrix3[(matrix3[0]<(binmax+1))&(matrix3[1]<(binmax+1))&(matrix3[2]<(binmax+1))]
	iterMax=3000
	if sys.argv[3]:
		z=np.zeros([binmax,binmax,binmax])
		for _,i in cis1.iterrows():
			z[i[0]-1,i[1]-1,i[2]-1]=i[3]
			z[i[0]-1,i[2]-1,i[1]-1]=i[3]    
			z[i[1]-1,i[0]-1,i[2]-1]=i[3]    
			z[i[1]-1,i[2]-1,i[0]-1]=i[3]    
			z[i[2]-1,i[0]-1,i[1]-1]=i[3]    
			z[i[2]-1,i[1]-1,i[0]-1]=i[3]    
		z=z.astype(float)
		old_bias = None
		m=z.shape[0]
		bias = np.ones((m,1,1))
		total_counts=z.sum()
		eps=0.0001
		for i in range(iterMax):
			sum_ts=z.sum(axis=0).sum(axis=0)
			dbias=sum_ts.reshape((m,1,1))
			dbias /= dbias[dbias != 0].mean()
			dbias**=(1/3)
			dbias[dbias == 0] = 1
			bias *= dbias
			z/=dbias
			z/=dbias.reshape((1,m,1))
			z/=dbias.reshape((1,1,m))
			bias *=((z.sum() / total_counts)**(1/3))
			z*=(total_counts/z.sum())
			if old_bias is not None and np.abs(old_bias - bias).sum() < eps:
				print(f'Iteration finishing at round ${i}')
				break
			old_bias=bias.copy()
		else:
			print(f'Iteration finishing at the MAX_INTERATIONS of ${iterMax}')
	out=[]
	for i in range(binmax):
		for j in range(i,binmax):
			for k in range(j,binmax):
				if z[i]!=0:
					out.append(''.join(map(str,[i,j,k,z[i,j,k]]))+'\n')
	with open(f'iced_${sys.argv[1]}','a') as o:
		for i in out:
			o.write(''.join(out))
	with open(f'iced_${sys.argv[1]}.bias','a') as o:
		for i,j in enumerate(bias):
			o.write('%s\t%s\n'%(i,j))
