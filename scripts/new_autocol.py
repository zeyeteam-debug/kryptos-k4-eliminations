"""Autokey substitution followed by columnar transposition, tested exactly against the K4 cribs.

Models (A = ST or KA alphabet, p_i = plaintext index, x_i = pre-transposition ciphertext index,
ct = columnar(x), i.e. ct[m[i]] = x_i with m[pt_index] = ct_index from gen.colperms):
  PT-auto   : x_i = p_i (+/-) p_{i-L}            (primer for i<L unknown -> no check there)
  CT-auto   : x_i = p_i (+/-) x_{i-L}            (classic ciphertext autokey, index before transposition)
  CT-post   : x_i = p_i (+/-) ct[i-L]            (key taken from the FINAL ciphertext index i-L)
  PT-refl   : x_i = p_i (+/-) p_{c-i}            (reflection about centre c, c in 40..160)
  CT-refl   : x_i = p_i (+/-) x_{c-i}
Rules: vig x=p+k, beau x=k-p, varb x=p-k.  A check is a crib position whose key letter is known;
a (perm,config) is a hit when every check matches and there are >=4 checks.
"""
import sys,time,random
import numpy as np
from stack import CT0,ANCH
from gen import colperms,ST,KA

RULES=('vig','beau','varb')
AN=np.array(sorted(ANCH))
N=97
LAGS=range(1,int(sys.argv[2]) if len(sys.argv)>2 else 31)
CENTRES=range(40,161)
MINCHK=4
WMAX=int(sys.argv[1]) if len(sys.argv)>1 else 7

def build_perms(wmax):
    Ms=[];labs=[]
    for W in range(2,wmax+1):
        M=colperms(N,W);Ms.append(M);h=len(M)//2
        labs+=[(W,k,'fwd') for k in range(h)]+[(W,k,'inv') for k in range(h)]
    return np.concatenate(Ms).astype(np.int64),labs

def combine(P,K,rule):
    if rule==0: return (P+K)%26
    if rule==1: return (K-P)%26
    return (P-K)%26

def run(ct,M,labs,verbose=False):
    """returns list of hits (family,alph,rule,param,nchecks,label) and expected false-positive count."""
    hits=[];expected=0.0
    cl=np.array([ord(c) for c in ct])
    for an,A in (('ST',ST),('KA',KA)):
        lut=np.zeros(128,dtype=np.int64);lut[[ord(c) for c in A]]=np.arange(26)
        v=lut[cl]                      # alphabet index of ct letter at each ct position
        CL=v[M]                        # CL[perm,i] = index of ct[m[i]]  (= x_i under that perm)
        Pidx={i:A.index(ANCH[i]) for i in ANCH}
        Pv=np.array([Pidx[i] for i in AN])
        for rule in range(3):
            # ---- PT-autokey and CT-autokey (pre / post) over lags ----
            for L in LAGS:
                # PT-auto: need i-L in cribs
                I=np.array([i for i in AN if i-L in ANCH]);
                if len(I)>=MINCHK:
                    X=combine(np.array([Pidx[i] for i in I]),np.array([Pidx[i-L] for i in I]),rule)
                    ok=np.all(CL[:,I]==X[None,:],axis=1)
                    expected+=len(M)*26.0**(-len(I))
                    for x in np.nonzero(ok)[0]: hits.append(('PT-auto',an,RULES[rule],L,len(I),labs[x]))
                # CT-auto (pre-transposition index): need i-L>=0
                I=AN[AN-L>=0]
                if len(I)>=MINCHK:
                    X=combine(np.array([Pidx[i] for i in I])[None,:],CL[:,I-L],rule)
                    ok=np.all(CL[:,I]==X,axis=1)
                    expected+=len(M)*26.0**(-len(I))
                    for x in np.nonzero(ok)[0]: hits.append(('CT-auto',an,RULES[rule],L,len(I),labs[x]))
                    # CT-post: key = final ciphertext letter at index i-L
                    X=combine(np.array([Pidx[i] for i in I]),v[I-L],rule)
                    ok=np.all(CL[:,I]==X[None,:],axis=1)
                    expected+=len(M)*26.0**(-len(I))
                    for x in np.nonzero(ok)[0]: hits.append(('CT-post',an,RULES[rule],L,len(I),labs[x]))
            # ---- reflection variants ----
            for c in CENTRES:
                I=np.array([i for i in AN if (c-i) in ANCH])
                if len(I)>=MINCHK:
                    X=combine(np.array([Pidx[i] for i in I]),np.array([Pidx[c-i] for i in I]),rule)
                    ok=np.all(CL[:,I]==X[None,:],axis=1)
                    expected+=len(M)*26.0**(-len(I))
                    for x in np.nonzero(ok)[0]: hits.append(('PT-refl',an,RULES[rule],c,len(I),labs[x]))
                I=AN[(c-AN>=0)&(c-AN<N)]
                if len(I)>=MINCHK:
                    X=combine(np.array([Pidx[i] for i in I])[None,:],CL[:,c-I],rule)
                    ok=np.all(CL[:,I]==X,axis=1)
                    expected+=len(M)*26.0**(-len(I))
                    for x in np.nonzero(ok)[0]: hits.append(('CT-refl',an,RULES[rule],c,len(I),labs[x]))
    return hits,expected

def plant_pt_auto(pt,A,L,primer,rule,W,order):
    x=[]
    for i,ch in enumerate(pt):
        k=A.index(primer[i]) if i<L else A.index(pt[i-L])
        x.append(int(combine(A.index(ch),k,rule)))
    return columnar(x,W,order,A)
def plant_ct_auto(pt,A,L,primer,rule,W,order):
    x=[]
    for i,ch in enumerate(pt):
        k=A.index(primer[i]) if i<L else x[i-L]
        x.append(int(combine(A.index(ch),k,rule)))
    return columnar(x,W,order,A)
def columnar(x,W,order,A):
    Lr=N//W;r=N%W;lens=[Lr+(1 if c<r else 0) for c in range(W)]
    idx=[]
    for c in order: idx+=[row*W+c for row in range(lens[c])]
    return ''.join(A[x[j]] for j in idx)

def summarize(tag,hits,expected):
    fam={}
    for h in hits: fam.setdefault(h[0],[]).append(h)
    best=sorted(hits,key=lambda h:-h[4])[:3]
    print('%s: hits=%d (expected null ~%.2f) by family %s; best %s'%(tag,len(hits),expected,{k:len(v) for k,v in fam.items()},best))

if __name__=='__main__':
    t0=time.time();M,labs=build_perms(WMAX);print('perms',len(M),'widths 2..%d (incl. inverses)'%WMAX)
    # positive control: PT-autokey KA L=7 vig, width 5 column order (3,0,4,1,2)
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    ctc=plant_pt_auto(pt,KA,7,'KRYPTOS',0,5,(3,0,4,1,2))
    h,e=run(ctc,M,labs);print('control PT-auto KA L=7 vig W5 order(3,0,4,1,2):',[x for x in h if x[0]=='PT-auto' and x[3]==7],'| other hits',len(h)-len([x for x in h if x[0]=='PT-auto' and x[3]==7]))
    ctc2=plant_ct_auto(pt,ST,4,'CLOK',1,6,(5,2,0,3,1,4))
    h,e=run(ctc2,M,labs);print('control CT-auto ST L=4 beau W6:',[x for x in h if x[0]=='CT-auto' and x[3]==4][:3],'| other hits',len(h)-len([x for x in h if x[0]=='CT-auto' and x[3]==4]))
    print('control time %.1fs'%(time.time()-t0))
    # K4
    t1=time.time();h,e=run(CT0,M,labs);summarize('K4',h,e);print('K4 run %.1fs'%(time.time()-t1))
    for x in sorted(h,key=lambda h:-h[4])[:12]: print('   ',x)
    # null: shuffled K4
    random.seed(7);NS=int(sys.argv[3]) if len(sys.argv)>3 else 20;tot=0;fam={};mx=0
    for s in range(NS):
        l=list(CT0);random.shuffle(l);hn,_=run(''.join(l),M,labs);tot+=len(hn)
        for x in hn: fam[x[0]]=fam.get(x[0],0)+1;mx=max(mx,x[4])
        if time.time()-t0>330: NS=s+1;break
    print('null: %d shuffles, mean hits/shuffle %.2f (expected %.2f), by family %s, max checks in a null hit %d'%(NS,tot/NS,e,fam,mx))
    print('total %.1fs'%(time.time()-t0))
