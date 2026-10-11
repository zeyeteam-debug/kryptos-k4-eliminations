"""A: substring-only systems (method changes part-way).  B: 4-layer PK-style stacks, 89 words.
usage: python3 new_subtext_4layer.py A|B|ctrlA|ctrlB"""
import numpy as np,itertools,random,sys,time
from stack import CT0,ANCH
from gen import run,colperms
from k3rot import rot
from layers3 import colidx,WORDS as W89base
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
T0=time.time()
def log(*a): print(f'[{time.time()-T0:6.0f}s]',*a,flush=True)
# ---------------- A ----------------
def perms_for(L):
    Ms=[np.arange(L,dtype=np.int16)[None,:]];labs=['id']
    for W in range(2,9):
        M=colperms(L,W);Ms.append(M);nf=len(M)//2
        labs+=[f'col{W}#{k}' for k in range(nf)]+[f'col{W}inv#{k}' for k in range(nf)]
    R=[];
    for w1 in range(2,L):
        for w2 in range(2,L):
            for d1 in (True,False):
                for d2 in (True,False):
                    idx=rot(L,w1,d1)[rot(L,w2,d2)]          # ct[n]=pt[idx[n]]
                    R.append(np.argsort(idx).astype(np.int16));labs.append(f'rot{w1}{"cw"[d1]}{w2}{"cw"[d2]}')
                    R.append(idx.astype(np.int16));labs.append(f'rot{w1}{"cw"[d1]}{w2}{"cw"[d2]}inv')
    Ms.append(np.array(R));return np.concatenate(Ms),labs
def runA(sub,anch,Ms,labs,CH=15000):
    h=[]
    for s in range(0,len(Ms),CH): h+=run(sub,anch,Ms[s:s+CH],labs[s:s+CH],Pmax=12)
    return h
def substrings():
    # run 1 (log new_subtext_A.log) did [0,74),[0,75),[0,82),[0,97); remaining set reduced for runtime
    # (dropped [14,74),[14,75),[14,82),[20,82),[21,82) which were not in the requested list)
    if len(sys.argv)>2 and sys.argv[2]=='rest':
        return [(14,97),(20,74),(20,75),(20,97),(21,74),(21,75),(21,97)]
    A=[0,14,20,21];B=[74,75,82,97];return [(a,b) for a in A for b in B]
def partA(ctrl=False):
    cache={}
    def P(L):
        if L not in cache: cache[L]=perms_for(L)
        return cache[L]
    if ctrl:
        rng=random.Random(5);a,b=21,97;L=b-a
        pt=list('THEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOESLANGLEYKNOWABOUTTHIS')[:L]
        anch={i-a:c for i,c in ANCH.items()}
        for i,c in anch.items(): pt[i]=c
        key='PALIMPSEST'[:8];cols=colperms(L,6)[123]          # m[pt]=ct
        sub=[None]*L
        for i,c in enumerate(pt): sub[cols[i]]=c
        sub=''.join(KA[(KA.index(c)+KA.index(key[j%8]))%26] for j,c in enumerate(sub))
        full=''.join(rng.choice(ST) for _ in range(a))+sub
        Ms,labs=P(L);h=runA(full[a:],anch,Ms,labs)
        log('ctrlA planted col6#123 + KA vig period 8 on [21,97):',len(h),'hits',[x for x in h if 'col6#123' in x[-1]][:3])
        return
    rng=random.Random(1);tot=0;totn=0
    for a,b in substrings():
        L=b-a;anch={i-a:c for i,c in ANCH.items()};Ms,labs=P(L)
        h=runA(CT0[a:b],anch,Ms,labs);s=list(CT0[a:b]);rng.shuffle(s)
        hn=runA(''.join(s),anch,Ms,labs);tot+=len(h);totn+=len(hn)
        log(f'[{a},{b}) L={L} perms={len(Ms)} K4 hits={len(h)} null hits={len(hn)}',h[:4],'| null:',hn[:2])
    log('A total K4',tot,'null',totn)
# ---------------- B ----------------
REQ=['KRYPTOS','PALIMPSEST','ABSCISSA','BERLIN','CLOCK','LAYERTWO','CLEPSYDRA','BEAMWORK','OCHRE','VERDIGRIS','UNDERLAY']
WB=[w for w in W89base if w not in ('ED','JIM','POOL','WATER','WOOD')]
WB=sorted(set(WB+REQ));assert len(WB)==89,len(WB)
def Topts(n,words):
    o=[('none',np.arange(n))]
    for w in words:
        if len(w)>1:
            idx=colidx(n,w);o.append(('T:'+w,np.argsort(idx)));o.append(('Ti:'+w,idx))   # pos[pt]=ct
    return o
SIGNS3=[(1,1,1),(1,1,-1),(1,-1,-1),(-1,-1,-1)]   # unordered multiset of signs suffices for sums of same-frame keys
def partB(ct,anch,words,AL,tag,want=None):
    n=len(ct);AN=np.array(sorted(anch));P=np.array([AL.index(anch[i]) for i in AN]);C=np.array([AL.index(c) for c in ct])
    V=[np.array([AL.index(c) for c in w]) for w in words];nW=len(words)
    def Kat(pos): return np.stack([v[pos%len(v)] for v in V])           # nW x 24
    TO=Topts(n,words);hits={'TQQQ':[],'QQQT':[],'TTQQ':[],'QTQT':[]}
    def tri(need,K):   # need = s1K1+s2K2+s3K3 (same frame); meet-in-middle
        out=[];d={}
        for i in range(nW):
            for j in range(i,nW):
                for sa,sb in ((1,1),(1,-1),(-1,1),(-1,-1)):
                    d.setdefault(((sa*K[i]+sb*K[j])%26).tobytes(),[]).append((i,sa,j,sb))
        for k in range(nW):
            for s in (1,-1):
                for (i,sa,j,sb) in d.get(((need-s*K[k])%26).tobytes(),[]): out.append((words[k],s,words[i],sa,words[j],sb))
        return out
    # T+Q+Q+Q : T on PT, then 3 Q in ct frame
    for tn,pos in TO:
        q=pos[AN];hits['TQQQ']+=[(tn,)+x for x in tri((C[q]-P)%26,Kat(q))]
    log(tag,'TQQQ done',len(hits['TQQQ']))
    # Q+Q+Q+T : Q's in PT frame then T
    Kp=Kat(AN)
    for tn,pos in TO:
        hits['QQQT']+=[(tn,)+x for x in tri((C[pos[AN]]-P)%26,Kp)]
    log(tag,'QQQT done',len(hits['QQQT']))
    # T+T+Q+Q and Q+T+Q+T
    for t1,p1 in TO:
        q1=p1[AN];K2=Kat(q1);dQT={}
        for j in range(nW):
            for s in (1,-1): dQT.setdefault(((s*K2[j])%26).tobytes(),[]).append((j,s))
        for t2,p2 in TO:
            q=p2[q1];need=(C[q]-P)%26
            # TTQQ: Q's at final ct positions q
            K=Kat(q);d={}
            for j in range(nW):
                for s in (1,-1): d.setdefault(((s*K[j])%26).tobytes(),[]).append((j,s))
            for s1 in (1,-1):
                L=(need[None,:]-s1*K)%26
                for i in range(nW):
                    for (j,s) in d.get(L[i].tobytes(),()): hits['TTQQ'].append((t1,t2,words[i],s1,words[j],s))
                # QTQT: Q1 on PT (Kp), T1, Q2 at q1 (K2), T2
                L=(need[None,:]-s1*Kp)%26
                for i in range(nW):
                    for (j,s) in dQT.get(L[i].tobytes(),()): hits['QTQT'].append((words[i],s1,t1,words[j],s,t2))
    log(tag,'TTQQ',len(hits['TTQQ']),'QTQT',len(hits['QTQT']))
    for k,v in hits.items():
        if want: v=[x for x in v if all(w in x for w in want)]
        log(tag,k,'hits',len(v),v[:4])
    return hits
def plant(order,AL,seed=3):
    """return ct with layers applied to a PT carrying the cribs"""
    rng=random.Random(seed)
    pt=list('ITWASTOTALLYINVISIBLEHOWSTHATPOSSIBLETHEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTED')[:97]
    for i,c in ANCH.items(): pt[i]=c
    x=[AL.index(c) for c in pt];n=97
    for kind,w,s in order:
        if kind=='Q': x=[(v+s*AL.index(w[i%len(w)]))%26 for i,v in enumerate(x)]
        else:
            idx=colidx(n,w)
            if s==-1: idx=np.argsort(idx)          # Ti
            x=[x[k] for k in idx]
    return ''.join(AL[v] for v in x)
if __name__=='__main__':
    m=sys.argv[1]
    if m=='A': partA()
    elif m=='ctrlA': partA(ctrl=True)
    elif m=='ctrlB':
        for AL,an in ((KA,'KA'),(ST,'ST')):
            for name,order,want in [
                ('TQQQ',[('T','PALIMPSEST',1),('Q','OCHRE',1),('Q','VERDIGRIS',-1),('Q','CLEPSYDRA',1)],['OCHRE','VERDIGRIS','CLEPSYDRA']),
                ('QQQT',[('Q','KRYPTOS',1),('Q','BEAMWORK',1),('Q','UNDERLAY',-1),('T','ABSCISSA',-1)],['KRYPTOS','BEAMWORK','UNDERLAY']),
                ('TTQQ',[('T','BERLIN',1),('T','LAYERTWO',-1),('Q','OCHRE',-1),('Q','CLOCK',1)],['OCHRE','CLOCK']),
                ('QTQT',[('Q','VERDIGRIS',1),('T','CLEPSYDRA',1),('Q','BEAMWORK',-1),('T','KRYPTOS',-1)],['VERDIGRIS','BEAMWORK'])]:
                ct=plant(order,AL);log('ctrlB',an,name,'planted',order)
                h=partB(ct,ANCH,WB,AL,f'ctrl-{an}-{name}',want=want)
    elif m=='B':
        rng=random.Random(7);s=list(CT0);rng.shuffle(s);NULL=''.join(s)
        for AL,an in ((KA,'KA'),(ST,'ST')):
            partB(CT0,ANCH,WB,AL,'K4-'+an)
            partB(NULL,ANCH,WB,AL,'NULL-'+an)
