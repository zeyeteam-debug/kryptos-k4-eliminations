import numpy as np,itertools,sys,time,random
from stack import CT0,ANCH
from spiral3 import spiral_order,perm_from_read,col_perm
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
WORDS=[w.strip().upper() for w in open('texts/dwyl/words.txt') if 2<=len(w.strip())<=30]
def colperms_all(N,W):
    out=[]
    for o in itertools.permutations(range(W)):
        idx=[];L=N//W;r=N%W
        for c in o: idx+=list(range(c,N,W))
        out.append(perm_from_read(idx,N))
    M=np.array(out);return np.concatenate([M,np.argsort(M,axis=1)])
def table(words,anch,AL):
    AN=np.array(sorted(anch));T={}
    for w in words:
        v=np.array([AL.index(ch) for ch in w])[AN%len(w)]
        for s in (1,-1): T.setdefault(((s*v)%26).astype(np.uint8).tobytes(),[]).append((w,s))
    return T
def run(ct,anch,routes,colWs,AL,TB):
    N=len(ct);AN=np.array(sorted(anch));P=np.array([AL.index(anch[i]) for i in AN]);C=np.array([AL.index(c) for c in ct])
    hits=[]
    cols=[(0,np.arange(N)[None,:])]+[(W,colperms_all(N,W)) for W in colWs]
    for rl,rm in routes:
        for W,Mc in cols:
            for order in ('R>T','T>R'):
                q=Mc[:,rm[AN]] if order=='R>T' else rm[Mc[:,AN]]
                D=((C[q]-P)%26).astype(np.uint8)
                for x in range(len(D)):
                    r=TB.get(D[x].tobytes())
                    if r: hits+=[(w,s,rl,W,x,order) for w,s in r]
    return hits
if __name__=='__main__':
    t=time.time()
    pk9=open('pk9.txt').read().replace('\n','').strip()
    p9='ISPENTTHEPASTMONTHWITHTHENEEDLEANDKNOTANDATLASTPELLEGRINSFINALMESSAGEHASBEENREVEALEDTOMEIWILLNOWSEALITFORYOUUNDEREVERYCIPHERIUSEDINTHISTESTAMENT'
    a9={i:p9[i] for i in list(range(21,34))+list(range(63,74))}
    TB9=table(WORDS,a9,KA);print('table built',len(TB9),round(time.time()-t),'s',flush=True)
    r9=[(('sp',12,'TR',1),perm_from_read(spiral_order(144,12,'TR',1),144))]
    h=run(pk9,a9,r9,[8],KA,TB9);print('PK9 control (370k-word Q, spiral, all width-8 orders):',[x for x in h if x[0]=='CLEPSYDRA'][:2],'total hits',len(h),round(time.time()-t),'s',flush=True)
    routes=[('none',np.arange(97))]
    for W in range(2,49):
        for corner in ('TL','TR','BL','BR'):
            for cw in (0,1):
                m=perm_from_read(spiral_order(97,W,corner,cw),97);routes+=[((W,corner,cw,'r'),m),((W,corner,cw,'w'),np.argsort(m))]
    colWs=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else [2,3,4,5,6]
    for an,AL in (('KA',KA),('ST',ST)):
        TB=table(WORDS,ANCH,AL)
        h=run(CT0,ANCH,routes,colWs,AL,TB);print(an,'K4: dictionary Q -> spiral/none -> columnar',colWs,': hits',len(h),h[:10],round(time.time()-t),'s',flush=True)
