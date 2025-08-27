#!/usr/bin/env python3
# 6G memory for initialing, 10G for 1.5M 3-way complexes. Using split
# command before using.
import sys
with open('/data/hic-pro-config/DpnII_hg38XX.bed') as f:
    d=[x.split() for x in f.readlines()]
with open('/data/hic-pro-config/hg38XX.chrom.sizes') as f:
	cs={y[0]:y[1] for x in f.readlines() for y in [x.split()]}
print('init...')
d=[[x[0],int(x[1]),int(x[2])+4,x[3]] for x in d]
d1={'%s_%s'%(x[0],x[1]):x[2] for x in d}
d2={'%s_%s'%(x[0],x[2]):x[1] for x in d}
d1n={'%s_%s'%(x[0],x[1]):x[3] for x in d}
d2n={'%s_%s'%(x[0],x[2]):x[3] for x in d}
d1k={x[0]:[] for x in d}
for x in d:
	d1k[x[0]].append([int(x[1]),int(x[2])])
da={x:{y//100000:[] for y in range(0,int(y),100000)} for x,y in cs.items()}
db={x:{y//100000:[] for y in range(0,int(y),100000)} for x,y in cs.items()}
del d
for i,j in d1k.items():
	for k in j:
		da[i][k[0]//100000].append(k[0])
		db[i][k[1]//100000].append(k[1])
for i,j in cs.items():
	for k in range(0,int(j)//100000):
		if len(da[i][k])==0:
			da[i][k].append(da[i][k-1][-1])
	for k in range(int(j)//100000,0,-1):
		if len(db[i][k])==0:
			db[i][k].append(db[i][k+1][0])
del cs
print('done')
def m3_to_bed(arg):
	file,arm,a,start=arg
	print(f'{arm}')
	with open('f_%s_%s.bed'%(file[:file.rindex('.')],arm),'w') as o1:
		with open('m_%s_%s.bed'%(file[:file.rindex('.')],arm),'w') as o2:
			for k,i in enumerate(a):
				if i[2] == '-':
					i1=int(i[1])
					bi=i1//100000
					while da[i[0]][bi][0] > i1:
						bi-=1
					for j in da[i[0]][bi][::-1]:
						if j<i1:
							p2=d1['%s_%s'%(i[0],j)]
							fn=d1n['%s_%s'%(i[0],j)]
							break
					else:
						j=da[i[0]][bi-1][-1]
						if j<i1:
							p2=d1['%s_%s'%(i[0],j)]
							fn=d1n['%s_%s'%(i[0],j)]
						else:
							print(i[0],j,i1)
							raise ValueError
					a1=[i[0],j,p2,k+start,fn,p2-j,i[2]]
				else:
					i2=int(i[1])
					bi=i2//100000
					while db[i[0]][bi][-1] < i2:
						bi+=1
					for j in db[i[0]][bi]:
						if j>i2:
							p1=d2['%s_%s'%(i[0],j)]
							fn=d2n['%s_%s'%(i[0],j)]
							break
					else:
						j=db[i[0]][bi+1][0]
						if j>i2:
							p1=d2['%s_%s'%(i[0],j)]
							fn=d2n['%s_%s'%(i[0],j)]
						else:
							print(i[0],j,i2)
							raise ValueError
					a1=[i[0],p1,j,k+start,fn,j-p1,i[2]]
				o1.write('\t'.join(map(str,a1))+'\n' )
				mi=(a1[1]+a1[2])//2
				b=[a1[0],mi,mi+1]+a1[3:]
				o2.write('\t'.join(map(str,b))+'\n')
	print(f'{arm} done')
if __name__=='__main__':
	from multiprocessing import Pool
	flag=int(sys.argv[1])
	num=0
	for file in sys.argv[2:]:
		if not flag:
			num=0
		print(f'{file}...')
		with open(file) as f:
			x=f.readline()
			w3=[[y.split(':')[1].split('|')] for y in x[:-1].split(';')]
			for x in f.readlines():
				for y,z in enumerate([y.split(':')[1].split('|') for y in x[:-1].split(';')]):
					w3[y].append(z)
		if '.' not in file:
			file+='.'
		with Pool(3) as pool:
			pool.map(m3_to_bed,[[file,x,w3[x],num] for x in range(len(w3))])
		num+=len(w3[0])
		print(f'{file} done')
