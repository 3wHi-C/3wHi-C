#!/usr/bin/env python3
import sys,gzip,os,re
'''def check(i):
    for j,x in enumerate(m.keys()):
        if sum([i[y]==x[y] for y in range(6)])==5:
            return j
    else:
        return -1'''
def c_seq(fq):
    c_fq={}
    for i,j in fq.items():
        if len(j)>16:
            if j.endswith('GATC'):
                c_fq.update({i:j})
            else:
                tmp=re.sub('GATC[GATC]{1,2}$|GATTC$|GAATC$|[ATC]ATC$|G[TCG]TC$|GA[GAC]C$|GAT[GAT]$|G[AT]C$|GAT$','GATC',j)
                if tmp.endswith('GATC'):
                    c_fq.update({i:tmp})
    return c_fq

def correct_seq(fq1,fq2):
    c_fq1=c_seq(fq1)
    c_fq2=c_seq(fq2)
    c_key=c_fq1.keys() & c_fq2.keys()
    c_fq1={x:c_fq1[x] for x in c_key}
    c_fq2={x:c_fq2[x] for x in c_key}
    return c_fq1,c_fq2

if __name__=='__main__':
    prefix=sys.argv[1]
    try:
        os.makedirs(f'hicpro/rawdata/{prefix}')
    except:
        pass
    with open(f'{prefix}.1.l.seq') as f:
        a=[x.split() for x in f.readlines()]
    u1=''
    u2=''
    fq1={}
    fq2={}
    m=dict(zip('ATCGN','TAGCN'))
    for i in a:
        b=''.join([m[x] for x in i[3][::-1]])
        fq1.update({i[0]:i[2][7:]})
        fq2.update({i[0]:b[7:]})
        u1+=i[0]+'\t'+i[2][:6]+'\n'
        u2+=i[0]+'\t'+b[:6]+'\n'
    fq1,fq2=correct_seq(fq1,fq2)
    with open(f'{prefix}_R1.umi','w') as f:
        f.write(u1)
    with open(f'{prefix}_R2.umi','w') as f:
        f.write(u2)
    with gzip.open(f'hicpro/rawdata/{prefix}/{prefix}_R1.fq.gz','wb') as f:
        f.write(''.join(['@%s\n%s\n+\n%s\n'%(x,y,'F'*len(y)) for x,y in fq1.items()]).encode())
    with gzip.open(f'hicpro/rawdata/{prefix}/{prefix}_R2.fq.gz','wb') as f:
        f.write(''.join(['@%s\n%s\n+\n%s\n'%(x,y,'F'*len(y)) for x,y in fq2.items()]).encode())
    '''with open(sys.argv[2]) as f:
        m={y[:-1]:x for x,y in enumerate(f.readlines())}
    with open(f'{prefix}_R1.umi') as f:
        u1={x[:x.index('\t')]:x[x.index('\t')+1:-1] for x in f.readlines()}
    with open(f'{prefix}_R2.umi') as f:
        u2={x[:x.index('\t')]:x[x.index('\t')+1:-1] for x in f.readlines()}
    with open(f'merged_{prefix}.fa') as f:
        fa={x[:x.index('\n')]:x[x.index('\n')+1:-1] for x in f.read().split('>')[1:]}
    with open(f'{prefix}.1.l.marker') as f:
        a=[x.split() for x in f.readlines()]
        a1={x[0]:x[3] for x in a if int(x[2])>21 and int(x[4])>21}
    with open(f'{prefix}.1.l.bin') as f:
        b=[x.split() for x in f.readlines()]
    bc={}
    for x in b:
        bc.update({x[0]:f'{x[0]}\t{x[1]}\t{int(x[2])-14}\t{int(x[3])-7} {int(x[4])>>7} {int(x[5])>>7} {int(x[6])>>7}\t{int(x[7])-7} {int(x[8])>>7} {int(x[9])>>7} {int(x[10])>>7}\n'})
    U1={x:u1[x] for x in a1.keys()}
    U2={x:u2[x] for x in a1.keys()}
    UA={y:m[x] if x in m.keys() else check(x) for y,x in U1.items()}
    UB={y:m[x] if x in m.keys() else check(x) for y,x in U2.items()}
    UA1={x:y for x,y in UA.items() if y!=-1}
    UB1={x:y for x,y in UB.items() if y!=-1}
    vu=(set(UA1.keys())&set(UB1.keys()))
    fq={y:{x:[] for x in range(64)} for y in range(64)}
    mV={x:y for x,y in a1.items() if x in vu}
    R1={'α':{x:[] for x in range(64)},'β':{x:[] for x in range(64)},'γ':{x:[] for x in range(64)}}
    R2={'α':{x:[] for x in range(64)},'β':{x:[] for x in range(64)},'γ':{x:[] for x in range(64)}}
    for x,y in mV.items():
        R1[y][UA1[x]].append(x)
        R2[y][UB1[x]].append(x)
        fq[UA1[x]][UB1[x]].append([x,fq1[x],fq2[x]])
    fil='adda\tα\t44\t1 0 0 0\t1 0 0 0\naddb\tβ\t48\t3 0 0 0\t3 0 0 0\naddg\tγ\t52\t5 0 0 0\t5 0 0 0\n'
    try:
        os.mkdir('match')
    except:
        pass
    os.chdir('match')
    for i in range(64):
        try:
            os.mkdir(str(i))
        except:
            pass
        os.chdir(str(i))
        with open('%s_a.1.l.bin'%i,'w') as o:
            o.write(fil+''.join([bc[x] for x in R1['α'][i] +R2['γ'][i] if x in bc.keys()]))
        with open('%s_b.1.l.bin'%i,'w') as o:
            o.write(fil+''.join([bc[x] for x in R1['β'][i] +R2['α'][i] if x in bc.keys()]))
        with open('%s_g.1.l.bin'%i,'w') as o:
            o.write(fil+''.join([bc[x] for x in R1['γ'][i] +R2['β'][i] if x in bc.keys()]))
        os.chdir('..')
    os.chdir('..')
    for i in range(64):
        for j in range(64):
            try:
                os.makedirs(f'hicpro_s/rawdata/{i}_{j}')
            except:
                pass
            with gzip.open(f'hicpro_s/rawdata/{i}_{j}/{i}_{j}_R1.fq.gz','wb') as f:
                f.write(''.join(['@%s\n%s\n+\n%s\n'%(x[0],x[1],'F'*len(x[1])) for x in fq[i][j]]).encode())
            with gzip.open(f'hicpro_s/rawdata/{i}_{j}/{i}_{j}_R2.fq.gz','wb') as f:
                f.write(''.join(['@%s\n%s\n+\n%s\n'%(x[0],x[2],'F'*len(x[2])) for x in fq[i][j]]).encode())'''
