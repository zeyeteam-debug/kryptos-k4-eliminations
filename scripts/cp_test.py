import sys,numpy as np
sys.path.insert(0,'.')
from stack import CT0,ANCH,KA,ST
from gen import run,colperms
# 1) Quagmire-II / Cyrillic-Projector form, exact: C = KA[(ST(P)+k)%26], k periodic. Also KA-plain/ST-cipher (Quag I).
def exact(ct,pa,ca,P):
    bad=0;m={}
    for i,p in ANCH.items():
        k=(ca.index(ct[i])-pa.index(p))%26
        if m.setdefault(i%P,k)!=k: bad+=1
    return bad
for nm,pa,ca in (("QII ST->KA (Cyrillic Projector form)",ST,KA),("QI KA->ST",KA,ST),("QIII KA",KA,KA),("Vig ST",ST,ST)):
    print(nm,[P for P in range(1,27) if exact(CT0,pa,ca,P)==0])
# keystream at cribs in CP form, as key letters in ST
ks={i:ST[(KA.index(CT0[i])-ST.index(p))%26] for i,p in sorted(ANCH.items())}
print('CP-form key letters (ST) at cribs:',''.join(ks[i] for i in sorted(ks)))
# MEDUSA-rule keys: word -> derive, test exactly
for w in ["MEDUSA","SHADOW","TEN","KRYPTOS","PALIMPSEST","ABSCISSA","BERLIN","CLOCK","ANTIPODES"]:
    for nm,pa,ca in (("QII",ST,KA),("QIII",KA,KA),("Vig",ST,ST)):
        hits=sum(ca[(pa.index(ANCH[i])+ST.index(w[i%len(w)]))%26]==CT0[i] for i in ANCH)
        hk=sum(ca[(pa.index(ANCH[i])+KA.index(w[i%len(w)]))%26]==CT0[i] for i in ANCH)
        if max(hits,hk)>=4: print('word',w,nm,hits,hk)
print('max crib hits from 9 words x3 forms printed only if >=4 (24 cribs; ~1 expected by chance)')
# 2) Transposition 2-8 + periodic <=12 via gen.run ('any' covers QII)
Ms=[];labs=[]
for W in range(2,9):
    M=colperms(97,W);Ms.append(M);labs+=[(W,j) for j in range(len(M))]
Ms=np.concatenate(Ms)
h=run(CT0,ANCH,Ms,labs,Pmax=12)
print('gen.run columnar 2-8 hits:',len(h),h[:10])
