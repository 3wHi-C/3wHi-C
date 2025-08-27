def func(a,b):
	f=1
	g=''
	for x in b:
		if x[0]-f<6:
			g+=x[2]+' '
		else:
			g+='%s %s '%(x[0]-f,x[2])
		f=x[1]
	if a-f<6:
		g=g[:-1]
	else:
		g+=str(a-f)
	return g
if __name__=='__main__':
	import sys
	prefix=sys.argv[1]
	with open('blast_%s/%s.blast'%(prefix,prefix)) as a:
		a1=a.readlines()
	a2=[x.split()[:-1] for x in a1]
	a3={x[0] for x in a2}
	a4={x:[] for x in a3}
	for x in a2:
		a4[x[0]].append([int(x[6]),int(x[7]),x[1][-1] if int(x[8])<int(x[9]) else x[1][-1].upper()])
	a5={x:sorted(y) for x,y in a4.items()}
	with open('%s.fa'%prefix) as b:
		b1=b.readlines()
	b2={b1[x][1:-1]:len(b1[x+1])-1 for x in range(0,len(b1),2)}
	a6={x:[] for x in b2.keys()}
	for x,y in a5.items():
		a6[x]=a5[x]
	c={}
	for x,y in b2.items():
		c.update({x:[y,func(y,a6[x])]})
	with open('%s.marker'%prefix,'w') as o:
		o.write(''.join('%s\t%s\t%s\n'%(x,y[0],y[1]) for x,y in c.items()))
	for i in [0,1,2,3]:
		c={}
		for x,y in b2.items():
			if len(a6[x])==i:
				c.update({x:[y,func(y,a6[x])]})
		with open('%s.%s.marker'%(prefix,i),'w') as o:
			o.write(''.join('%s\t%s\t%s\n'%(x,y[0],y[1]) for x,y in c.items()))
		if i==1:
			c={x:y for x,y in c.items() if len(y[1].split())==3}
			with open('%s.%s.l.marker'%(prefix,i),'w') as o:
				o.write(''.join('%s\t%s\t%s\n'%(x,y[0],y[1]) for x,y in c.items()))
