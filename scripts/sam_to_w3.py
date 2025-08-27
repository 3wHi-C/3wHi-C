#!/usr/bin/env python3
import sys
from collections import Counter
from random import sample
from matplotlib import pyplot as plt
from functools import reduce
import re
def check(i):
    for j,x in enumerate(m.keys()):
        if sum([i[y]==x[y] for y in range(6)])==5:
            return j
    else:
        return -1
def makeBias(pos,bias):
    a=pos.split('|')
    return ''.join([a[0],str(int(a[1])+bias if a[2]=='+' else int(a[1])-bias),a[2]])
def correct_len(i):
    if len(i)<5:
        return int(i[:-1])
    else:
        ad=re.findall('\d+(?=[MD])',i)
        de=re.findall('\d+(?=I)',i)
        return sum([int(x) for x in ad])-sum([int(x) for x in de])
def get_w3(UA1,UB1,marker,sp,piece_len=0):
    if piece_len:
        sp1=sample(range(len(sp)//2),piece_len//2)
        sp=[sp[y] for x in sp1 for y in [x*2,x*2+1]]
        del sp1
    #else:
        #piece_len=len(sp)
    sam={}
    for x in range(0,len(sp),2):
        flag1=int(sp[x][1])>>4&1
        flag2=int(sp[x+1][1])>>4&1
        len1=round((correct_len(sp[x][5])-4)/5)*5+3 if flag1 else -(correct_len(sp[x][5])+2)%5
        len2=round((correct_len(sp[x+1][5])-4)/5)*5+3 if flag2 else -(correct_len(sp[x][5])+2)%5
        sam.update({sp[x][0]:[[sp[x][2],int(sp[x][3])+len1,'-' if flag1 else '+'],[sp[x+1][2],int(sp[x+1][3])+len2,'-' if flag2 else '+']]})
    
    vu=(set(UA1.keys())&set(UB1.keys()))
    sI=set(sam)&vu
    SS={f'{UA1[x]};{UB1[x]};{sam[x][0][0]}|{sam[x][0][1]}|{sam[x][0][2]};{sam[x][1][0]}|{sam[x][1][1]}|{sam[x][1][2]}':x for x in sI}
    Ss={y:[int(z[0]),int(z[1]),z[2],z[3]] for x,y in SS.items() for z in [x.split(';')]}
    lSs=len(Ss)
    sV={x:y for x,y in marker.items() if x in sI}
    Rs1={'α':{x:[] for x in range(64)},'β':{x:[] for x in range(64)},'γ':{x:[] for x in range(64)}}
    Rs2={'α':{x:[] for x in range(64)},'β':{x:[] for x in range(64)},'γ':{x:[] for x in range(64)}}
    for x,y in Ss.items():
        Rs1[sV[x]][y[0]].append([y[2],f'{y[1]}:{y[3]}'])
        Rs2[sV[x]][y[1]].append([y[3],f'{y[0]}:{y[2]}'])
    RS1={'α':{x:{} for x in range(64)},'β':{x:{} for x in range(64)},'γ':{x:{} for x in range(64)}}
    RS2={'α':{x:{} for x in range(64)},'β':{x:{} for x in range(64)},'γ':{x:{} for x in range(64)}}
    singleDupCount='\tR1\tR2\n'
    LG1=[]
    Lg1=[]
    LG2=[]
    Lg2=[]
    for i in 'αβγ':
        for j in range(64):
            g=Rs1[i][j]
            G={i[0]:i[1] for i in g}
            G1={x for x,y in Counter([x[0] for x in g]).items() if y==1}
            RS1[i][j]={x:G[x] for x in G1}
            lg=len(g)
            lG=len(G1)
            Lg1.append(lg)
            LG1.append(lG)
            if lg==0:
                ra='NA'
            else:
                ra=lG/lg
            singleDupCount+=f'{i}_{j}\t{lG}/{lg}({ra})'
            g=Rs2[i][j]
            G={x[0]:x[1] for x in g}
            G1={x for x,y in Counter([x[0] for x in g]).items() if y==1}
            RS2[i][j]={x:G[x] for x in G1}
            lg=len(g)
            lG=len(G1)
            Lg2.append(lg)
            LG2.append(lG)
            if lg==0:
                ra='NA'
            else:
                ra=lG/lg
            singleDupCount+=f'\t{lG}/{lg}({ra})\n'
    with open(f'{piece_len}_single_duplicate.count','w') as o:
        if sum(Lg1)==0:
            ra1='NA'
        else:
            ra1=sum(LG1)/sum(Lg1)
        if sum(Lg2)==0:
            ra2='NA'
        else:
            ra2=sum(LG2)/sum(Lg2)
        o.write(singleDupCount+f'total\t{sum(LG1)}/{sum(Lg1)}({ra1})\t{sum(LG2)}/{sum(Lg2)}({ra2})\n')
    with open('A_valid.m2','w') as o:
        o.write(''.join([f'{x}:{z};{z1}\n' for x,y in RS1['α'].items() for z,z1 in y.items()]))
    with open('B_valid.m2','w') as o:
        o.write(''.join([f'{x}:{z};{z1}\n' for x,y in RS1['β'].items() for z,z1 in y.items()]))
    with open('G_valid.m2','w') as o:
        o.write(''.join([f'{x}:{z};{z1}\n' for x,y in RS1['γ'].items() for z,z1 in y.items()]))
    del Rs1,Rs2,lg,LG1,LG2,Lg1,Lg2
    mm={'α':{x:{} for x in range(64)},'β':{x:{} for x in range(64)},'γ':{x:{} for x in range(64)}}
    for i in range(64):
        Zs1=RS1['β'][i]
        Zs2=RS2['α'][i]
        ma=set(Zs1.keys()) & set(Zs2.keys())
        mm['α'][i]={f'{i}:{x};{Zs2[x]};{Zs1[x]}' for x in ma}
        Zs1=RS1['γ'][i]
        Zs2=RS2['β'][i]
        ma=set(Zs1.keys()) & set(Zs2.keys())
        mm['β'][i]={f'{Zs2[x]};{Zs1[x]};{i}:{x}' for x in ma}
        Zs1=RS1['α'][i]
        Zs2=RS2['γ'][i]
        ma=set(Zs1.keys()) & set(Zs2.keys())
        mm['γ'][i]={f'{Zs1[x]};{i}:{x};{Zs2[x]}' for x in ma}
    del RS1,RS2,Zs1,Zs2,ma
    AA={y for x in mm['α'].values() for y in x}
    BB={y for x in mm['β'].values() for y in x}
    GG={y for x in mm['γ'].values() for y in x}
    del mm
    W3=(AA&BB)|(BB&GG)|(GG&AA)
    A=AA-W3
    B=BB-W3
    G=GG-W3
    with open(f'A_{piece_len}.m3','w') as o:
        o.write('\n'.join(A)+'\n')
    with open(f'B_{piece_len}.m3','w') as o:
        o.write('\n'.join(B)+'\n')
    with open(f'G_{piece_len}.m3','w') as o:
        o.write('\n'.join(G)+'\n')
    with open(f'W3_{piece_len}.m3','w') as o:
        o.write('\n'.join(W3)+'\n')
    func=lambda x,y,z:[x[y],x[z]]
    RA=[func(x.split(';'),1,2) for x in A]
    RB=[func(x.split(';'),0,1) for x in B]
    RG=[func(x.split(';'),2,0) for x in G]
    RM=RA+RB+RG
    with open(f'supposed_{piece_len}.m2','w') as o:
        o.write(''.join([';'.join(x)+'\n' for x in RM]))
    return f'{piece_len}\t{lSs}\t{len(A)+len(B)+len(G)+len(W3)}\t{len(W3)}\t{len(AA&BB)}\t{len(BB&GG)}\t{len(GG&AA)}\t{len(AA&BB&GG)}\n'


if __name__=='__main__':
    prefix=sys.argv[1]
    with open('%s.1.l.marker'%prefix) as f:
        a=[x.split() for x in f.readlines()]
        a1={x[0]:x[3] for x in a if int(x[2])>21 and int(x[4])>21}
    del a
    with open('../UMI.map') as f:
        m={y[:-1]:x for x,y in enumerate(f.readlines())}
    with open('%s_R1.umi'%prefix) as f:
        u1={x[:x.index('\t')]:x[x.index('\t')+1:-1] for x in f.readlines()}
    with open('%s_R2.umi'%prefix) as f:
        u2={x[:x.index('\t')]:x[x.index('\t')+1:-1] for x in f.readlines()}
    U1={x:u1[x] for x in a1.keys()}
    U2={x:u2[x] for x in a1.keys()}
    UA={y:m[x] if x in m.keys() else check(x) for y,x in U1.items()}
    UB={y:m[x] if x  in m.keys() else check(x) for y,x in U2.items()}
    UA1={x:y for x,y in UA.items() if y!=-1}
    UB1={x:y for x,y in UB.items() if y!=-1}
    del m,UA,UB
    with open('%s_hg38XX.bwt2pairs_interaction.sam'%prefix) as f:
        sp=[x.split() for x in f.readlines()]
    
    num=[]
    if len(sys.argv)==3:
        for i in range(500000,len(sp),500000):
            num.append(get_w3(UA1,UB1,a1,sp,i))
    num.append(get_w3(UA1,UB1,a1,sp))
    with open('match_count.txt','w') as o:
        o.write(''.join(num))
    

