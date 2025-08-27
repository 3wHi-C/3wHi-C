#!/usr/bin/env python3
if __name__=='__main__':
	import sys,os
	prefix=sys.argv[1]
	f=lambda x,y:[x[0],int(x[1]),x[6],x[4],int(x[5]),y]
	z=0
	with open('%s_0.bed'%prefix) as f0:
		with open('%s_1.bed'%prefix) as f1:    
			with open('%s.validpairs'%prefix,'w') as o:
				if os.path.exists('%s_2.bed'%prefix):
					with open('%s_2.bed'%prefix) as f2:
						l0=f0.readline()
						l1=f1.readline()
						l2=f2.readline()
						while l0:
							w0=l0.split()
							w1=l1.split()
							w2=l2.split()
							l=sorted([f(w0,0),f(w1,1),f(w2,2)])
							o.write('%s.%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t42\t42\n'\
									%(z,l[0][-1],l[1][-1],l[0][0],l[0][1],l[0][2],l[1][0],l[1][1],l[1][2],\
									 l[0][4]+l[1][4],l[0][3],l[1][3]))
							o.write('%s.%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t42\t42\n'\
									%(z,l[0][-1],l[2][-1],l[0][0],l[0][1],l[0][2],l[2][0],l[2][1],l[2][2],\
									 l[0][4]+l[2][4],l[0][3],l[2][3]))
							o.write('%s.%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t42\t42\n'\
									%(z,l[1][-1],l[2][-1],l[1][0],l[1][1],l[1][2],l[2][0],l[2][1],l[2][2],\
									 l[1][4]+l[2][4],l[1][3],l[2][3]))
							l0=f0.readline()
							l1=f1.readline()                    
							l2=f2.readline()                    
							z+=1
				else:
					l0=f0.readline()
					l1=f1.readline()
					while l0:
						w0=l0.split()
						w1=l1.split()
						l=sorted([f(w0,0),f(w1,1)])
						o.write('%s.%s%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t42\t42\n'\
								%(z,l[0][-1],l[1][-1],l[0][0],l[0][1],l[0][2],l[1][0],l[1][1],l[1][2],\
									 l[0][4]+l[1][4],l[0][3],l[1][3]))
						l0=f0.readline()
						l1=f1.readline() 
						z+=1
