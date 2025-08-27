#!/usr/bin/env python3
from math import log2
class bin_fa:
	def __init__(self):
		self.map1={'A':0,'T':0,'C':1,'G':1,'N':1}
		self.map2={'A':0,'T':1,'C':0,'G':1,'N':1}
		self.map2r={'A':1,'T':0,'C':1,'G':0,'N':0}
		self.N={'A':1,'T':1,'C':1,'G':1,'N':0}
	
	def translate_PE(self,seq1,seq2):
		if self._check(seq1[-4:]) and self._check(seq2[:4]) and len(seq1)>12 and len(seq2)>12:
			return self.translate_r(seq1[:-4]+'GATC'),self.translate('GATC'+seq2[4:])
		return 0
		
	def _check(self,seq):
		a=self._translate(seq)
		b=(a[0]^9)|(a[1]^10)
		if b and not log2(b).is_integer():
			return 0
		return 1
	
	def _translate(self,seq):
		bin1=0
		bin2=0
		for i in seq:
			bin1<<=1
			bin1+=self.map1[i]
			bin2<<=1
			bin2+=self.map2[i]
		return (bin1,bin2)
	def translate(self,seq):
		bin1=0
		bin2=0
		N0=0
		f=0
		for i in seq:
			bin1<<=1
			bin1+=self.map1[i]
			bin2<<=1
			bin2+=self.map2[i]
			N0<<=1
			N0+=self.N[i]
		return (len(seq),bin1,bin2,N0)
	def translate_r(self,seq):
		bin1=0
		bin2=0
		N0=0
		f=0
		for i in seq[::-1]:
			bin1<<=1
			bin1+=self.map1[i]
			bin2<<=1
			bin2+=self.map2r[i]
			N0<<=1
			N0+=self.N[i]
			f+=1
			f%=5
		return (len(seq),bin1,bin2,N0)
		
def func(input_dict):
	from itertools import islice
	fa=bin_fa()
	in_keys=islice(input_dict[0].keys(),int(len(input_dict[0])*input_dict[1]/input_dict[2]),int(len(input_dict[0])*(input_dict[1]+1)/input_dict[2]))
	r={}
	for i in in_keys:
		result=fa.translate_PE(input_dict[0][i][1],input_dict[0][i][2])
		if result:
			r.update({i:[input_dict[0][i][0],result]})
	return r
if __name__ == '__main__':
	from multiprocessing import Pool
	from functools import reduce
	import sys
	prefix=sys.argv[1]
	if len(sys.argv)==2:
		t=10
	else:
		t=int(sys.argv[2])
	with open('%s.1.l.marker'%prefix) as a:
		a1=[x.split() for x in a.readlines()]
	a2={x[0]:x[1:] for x in a1}
	del a1
	with open('merged_%s.fa'%prefix) as b:
		b1=b.readlines()
	b2={b1[x][1:-1]:b1[x+1][:-1] for x in range(0,len(b1),2)}
	c={}
	for i,j in a2.items():
		c.update({i:[j[0],b2[i][:int(j[1])],b2[i][-int(j[3]):]]})
	del b2
	with open('%s.1.l.seq'%prefix,'w') as o:
		o.write(''.join(['%s\t%s\n'%(x,'\t'.join(y)) for x,y in c.items()]))
	with Pool(t) as pool:
		c=pool.map(func,[[c,x,t] for x in range(t)])
	d={}
	for x in c:
		d.update(x)
	del c
	with open('%s.1.l.bin'%prefix,'w') as o:
		o.write(''.join(['%s\t%s\t%s\t%s\t%s\n'%(x,a2[x][2],y[0],
		 ' '.join([str(x) for x in y[1][0]]),' '.join([str(x) for x in y[1][1]]))
		 for x,y in d.items()]))
