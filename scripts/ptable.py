import numpy as np,random,math
from stack import CT0,ANCH
from gen import run
from stackT import scan as stackscan
from spiral3 import perm_from_read
exec(open('board.py').read().split('# K4 cell')[0])
from sculpt import TAB
KA=KRY;ST=STD
def layout(kind):
    pos={}
    if kind=='18':     # standard: f-block in rows 9,10 cols 3..17 (La..Lu, Ac..Lr)
        pos[1]=(1,1);pos[2]=(1,18)
        for p,z0 in ((2,3),(3,11)):
            pos[z0]=(p,1);pos[z0+1]=(p,2)
            for k in range(6): pos[z0+2+k]=(p,13+k)
        for p,z0 in ((4,19),(5,37)):
            for k in range(18): pos[z0+k]=(p,1+k)
        for p,z0,f in ((6,55,9),(7,87,10)):
            pos[z0]=(p,1);pos[z0+1]=(p,2)
            for k in range(15): pos[z0+2+k]=(f,3+k)
            for k in range(15): pos[z0+17+k]=(p,4+k)
    else:              # 32-column long form
        pos[1]=(1,1);pos[2]=(1,32)
        for p,z0 in ((2,3),(3,11)):
            pos[z0]=(p,1);pos[z0+1]=(p,2)
            for k in range(6): pos[z0+2+k]=(p,27+k)
        for p,z0 in ((4,19),(5,37)):
            pos[z0]=(p,1);pos[z0+1]=(p,2)
            for k in range(16): pos[z0+2+k]=(p,17+k)
        for p,z0 in ((6,55),(7,87)):
            for k in range(32): pos[z0+k]=(p,1+k)
    return pos
def routes():
    out=[]
    for kind in ('18','32'):
        pos=layout(kind)
        for s in range(1,119-96):
            Z=list(range(s,s+97))
            for rd in ('rows','cols','rows_rev','cols_rev'):
                key={'rows':lambda z:pos[z],'cols':lambda z:(pos[z][1],pos[z][0]),'rows_rev':lambda z:(-pos[z][0],-pos[z][1]),'cols_rev':lambda z:(-pos[z][1],-pos[z][0])}[rd]
                order=sorted(range(97),key=lambda k:key(Z[k]))     # letter k sits on element Z[k]; read grid
                m=perm_from_read(order,97)
                out.append(((kind,s,rd,'r'),m));out.append(((kind,s,rd,'w'),np.argsort(m)))
    return out
def keystreams():
    for kind in ('18','32'):
        pos=layout(kind)
        for gname,G in (('tableau',TAB),('panel',LINES)):
            R=len(G);Wd=min(len(r) for r in G)
            for r0 in range(R):
                for c0 in range(Wd):
                    for s in range(1,119-96):
                        ks=[]
                        for z in range(s,s+97):
                            r,c=pos[z];rr,cc=r0+r-1,c0+c-1
                            if rr>=R or cc>=len(G[rr]) or not G[rr][cc].isalpha(): ks=None;break
                            ks.append(G[rr][cc])
                        if ks: yield (kind,gname,r0,c0,s),''.join(ks)
if __name__=='__main__':
    rt=routes();Ms=np.array([m.astype(np.int16) for _,m in rt]);labs=[l for l,_ in rt];print('periodic-table routes',len(Ms))
    print('route + any periodic key:',run(CT0,ANCH,Ms,labs,Pmax=24,RTHR=5,ANYSTRONG=3))
    print('route + two stacked keys:',stackscan(CT0,Ms,labs,amax=10,bmax=10))
    from dictsearch import table,run as drun,WORDS
    for an,AL in (('KA',KA),('ST',ST)):
        TB=table(WORDS,ANCH,AL);print(an,'dictionary Q + periodic-table route (+columnar 2-5):',drun(CT0,ANCH,[(l,m) for l,m in rt],[2,3,4,5],AL,TB)[:5])
    AN=sorted(ANCH);best=0;hits=[];n=0
    for lab,ks in keystreams():
        for an,A in (('KA',KA),('ST',ST)):
            for rule in range(3):
                m=0
                for i in AN:
                    k=A.index(ks[i]);p=A.index(ANCH[i]);m+=A[[(p+k)%26,(k-p)%26,(p-k)%26][rule]]==CT0[i]
                n+=1;best=max(best,m)
                if m==24: hits.append((lab,an,rule))
    print('periodic-table-shaped mask over tableau/panel as KEY: tested',n,'best',best,'/24 hits',hits)
