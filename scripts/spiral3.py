import numpy as np,json,sys,time
from stack import CT0,ANCH
from biglist import WORDS
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
def spiral_order(N,W,corner,cw):
    R=(N+W-1)//W
    starts={'TL':(0,0),'TR':(0,W-1),'BL':(R-1,0),'BR':(R-1,W-1)}
    seqs={('TL',1):'RDLU',('TL',0):'DRUL',('TR',1):'DLUR',('TR',0):'LDRU',('BL',1):'URDL',('BL',0):'RULD',('BR',1):'LURD',('BR',0):'ULDR'}
    D={'D':(1,0),'L':(0,-1),'U':(-1,0),'R':(0,1)};dirs=seqs[(corner,cw)]
    r,c=starts[corner];seen=set();out=[];d=0
    for _ in range(R*W):
        seen.add((r,c))
        if r*W+c<N: out.append(r*W+c)
        for t in range(4):
            dr,dc=D[dirs[d%4]];nr,nc=r+dr,c+dc
            if 0<=nr<R and 0<=nc<W and (nr,nc) not in seen: break
            d+=1
        r,c=nr,nc
    return out
def perm_from_read(order,N):      # out[k]=in[order[k]] ; return m[in]=k
    m=np.empty(N,dtype=np.int32);m[np.array(order)]=np.arange(N);return m
def order_of(k): return sorted(range(len(k)),key=lambda i:(k[i],i))
def col_perm(N,key):
    idx=[];w=len(key)
    for c in order_of(key): idx+=list(range(c,N,w))
    return perm_from_read(idx,N)
def build(N,widths,words):
    S=[]
    for W in widths:
        for corner in ('TL','TR','BL','BR'):
            for cw in (0,1):
                m=perm_from_read(spiral_order(N,W,corner,cw),N)
                S.append(((W,corner,cw,'read'),m));S.append(((W,corner,cw,'write'),np.argsort(m)))
    T=[('none',np.arange(N))]
    for w in words:
        if len(w)>1:
            m=col_perm(N,w);T.append(('T:'+w,m));T.append(('Tinv:'+w,np.argsort(m)))
    return S,T
def search(ct,anch,S,T,words,AL):
    N=len(ct);AN=np.array(sorted(anch));P=np.array([AL.index(anch[i]) for i in AN]);C=np.array([AL.index(c) for c in ct])
    table={}
    for w in words:
        v=np.array([AL.index(ch) for ch in w])[AN%len(w)]
        for s in (1,-1): table.setdefault(((s*v)%26).astype(np.int64).tobytes(),[]).append((w,s))
    hits=[]
    for sl,sm in S:
        for tl,tm in T:
            for order in ('S>T','T>S'):
                q=tm[sm[AN]] if order=='S>T' else sm[tm[AN]]
                d=((C[q]-P)%26).astype(np.int64).tobytes()
                for w,s in table.get(d,[]): hits.append((w,s,sl,tl,order))
    return hits
if __name__=='__main__':
    
    pk9=open('pk9.txt').read().replace('\n','').strip();p9='ISPENTTHEPASTMONTHWITHTHENEEDLEANDKNOTANDATLASTPELLEGRINSFINALMESSAGEHASBEENREVEALEDTOMEIWILLNOWSEALITFORYOUUNDEREVERYCIPHERIUSEDINTHISTESTAMENT';N9=len(pk9)
    a9={i:p9[i] for i in list(range(21,34))+list(range(63,74))}
    words=WORDS+['CLEPSYDRA','BEAMWORK']
    S,T=build(N9,[12],words)
    t=time.time();h=search(pk9,a9,S,T,words,KA);print('PK9 control:',h[:3],len(h),round(time.time()-t),'s',flush=True)
    S,T=build(97,range(2,49),WORDS);print('K4 spiral variants',len(S),'transpositions',len(T),flush=True)
    for an,AL in (('KA',KA),('ST',ST)):
        t=time.time();h=search(CT0,ANCH,S,T,WORDS,AL);print(an,'K4 Q->spiral->T hits:',len(h),h[:5],round(time.time()-t),'s',flush=True)
