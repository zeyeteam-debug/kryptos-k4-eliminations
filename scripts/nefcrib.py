import numpy as np,random,string
from stack import CT0,ANCH
from gen import run
W=['NEFERTITI','NOFRETETE','QUEENNEFERTITI','BUSTOFNEFERTITI','AKHENATEN','AMARNA','THUTMOSE','BORCHARDT','EGYPT','EGYPTIAN','CAIRO','LUXOR','TUTANKHAMUN','CARTER']
AN0=set(ANCH)
def slots(L): return [s for s in range(0,98-L) if not (set(range(s,s+L))&AN0)]
def test(w):
    n=0;key=0;anyp=[]
    for s in slots(len(w)):
        an=dict(ANCH);an.update({s+j:c for j,c in enumerate(w)});n+=1
        H=run(CT0,an,np.arange(97,dtype=np.int16)[None,:],[(w,s)],Pmax=48,RTHR=4,ANYSTRONG=3)
        key+=sum(1 for h in H if h[0]=='key');anyp+=[(h[2],s) for h in H if h[0]=='any']
    return n,key,anyp
rng=random.Random(5)
for w in W:
    n,k,a=test(w)
    rn=0;rk=0;ra=0
    for _ in range(20):
        r=''.join(rng.choice(string.ascii_uppercase) for _ in w);x=test(r);rn+=x[0];rk+=x[1];ra+=len(x[2])
    print(f'{w:16s} positions {n:3d} | ST/KA periodic<=48 consistent: {k} | any-alphabet periodic consistent: {len(a)} {a[:4]} | random same-length words: per-position rate {ra/max(rn,1):.3f}')
