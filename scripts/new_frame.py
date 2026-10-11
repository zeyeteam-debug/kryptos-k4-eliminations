"""Test the 'indicator' hypothesis: OBKR (first 4) / KCAR (last 4) are not message,
or a '?' cell precedes OBKR. Columnar transposition widths 2..W + periodic key via gen.run."""
import numpy as np,random,time,sys
from stack import CT0,ANCH
from gen import run,colperms
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
Pmax=12;BATCH=20000
WMAX=int(sys.argv[1]) if len(sys.argv)>1 else 8

def shifted(d):
    return {i+d:c for i,c in ANCH.items()}
VAR={
 'a_drop_first4': (CT0[4:],       shifted(-4)),
 'b_drop_last4':  (CT0[:93],      shifted(0)),
 'c_drop_both':   (CT0[4:93],     shifted(-4)),
 'd_prepend_q':   ('?'+CT0,       shifted(+1)),
 'e_prepend_q_drop_last4': ('?'+CT0[:93], shifted(+1)),
}
def null_of(ct,seed):
    idx=[i for i,c in enumerate(ct) if c!='?'];l=[ct[i] for i in idx]
    random.Random(seed).shuffle(l);out=list(ct)
    for i,c in zip(idx,l): out[i]=c
    return ''.join(out)
def sweep(ct,anch,Ws,tag):
    N=len(ct);H=[];nperm=0;t=time.time()
    for W in Ws:
        Ms=colperms(N,W);nperm+=len(Ms)
        for c0 in range(0,len(Ms),BATCH):
            L=[(W,i) for i in range(c0,min(c0+BATCH,len(Ms)))]
            H+=run(ct,anch,Ms[c0:c0+BATCH],L,Pmax=Pmax)
    print(f'  {tag}: N={N} widths {Ws[0]}..{Ws[-1]} perms={nperm} hits={len(H)} {H[:4]} ({time.time()-t:.1f}s)',flush=True)
    return H,nperm

Ws=list(range(2,WMAX+1))
# ---------- positive control on variant (a) frame ----------
ct_a,an_a=VAR['a_drop_first4'];N=len(ct_a)
rng=random.Random(2024);pt=[KA[rng.randrange(26)] for _ in range(N)]
for i,c in an_a.items(): pt[i]=c
sub=[KA[(KA.index(c)+KA.index('LUCID'[i%5]))%26] for i,c in enumerate(pt)]   # KA Vigenere, key by pt position
M6=colperms(N,6);k=37;m=M6[k];ctrl=['?']*N
for i in range(N): ctrl[m[i]]=sub[i]
ctrl=''.join(ctrl)
h=run(ctrl,an_a,M6,[(6,i) for i in range(len(M6))],Pmax=Pmax)
key_hits=[x for x in h if x[0]=='key']
print('control (plant: W6 perm #%d, KA vig LUCID p5): key hits'%k,key_hits[:3],'| any-alpha hits',len(h)-len(key_hits),flush=True)
assert any(x[1]=='KA' and x[2]==0 and x[3]=='pt' and x[4]==5 and x[6]==(6,k) for x in key_hits),'control NOT found'
Hc,_=sweep(ctrl,an_a,Ws,'control full sweep')
print('  control full-sweep key hits:',[x for x in Hc if x[0]=='key'],flush=True)

# ---------- K4 variants + shuffled nulls ----------
tot=0
for name,(ct,anch) in VAR.items():
    print(name,flush=True)
    Hk,n=sweep(ct,anch,Ws,'K4  ');tot+=n
    Hn,_=sweep(null_of(ct,11),anch,Ws,'null')
    print(f'  => K4 hits {len(Hk)} (key {sum(x[0]=="key" for x in Hk)}, any {sum(x[0]=="any" for x in Hk)}) | null hits {len(Hn)}',flush=True)
print('total K4 perm-tests',tot)
