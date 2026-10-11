import sys,random,itertools,numpy as np
sys.path.insert(0,'.')
from stack import CT0,ANCH,KA,ST
from gen import colperms
AN=np.array(sorted(ANCH));PR=np.array(list(itertools.combinations(range(len(AN)),2)))
def scan(ct,M):
    out=[]
    cl=np.array([[ST,KA].__len__() and 0]*0)
    for pa_n,pa,ca_n,ca in (("ST",ST,"KA",KA),("KA",KA,"ST",ST),("ST",ST,"ST",ST),("KA",KA,"KA",KA)):
        C=np.array([ca.index(c) for c in ct])[M[:,AN]];P=np.array([pa.index(ANCH[i]) for i in AN])
        for s in (1,-1):
            K=(C-s*P)%26;same=K[:,PR[:,0]]==K[:,PR[:,1]]
            for Per in range(1,13):
                for kind,R in (('pt',np.broadcast_to(AN%Per,M[:,AN].shape)),('ct',M[:,AN]%Per)):
                    need=R[:,PR[:,0]]==R[:,PR[:,1]]
                    ok=np.all(same|~need,axis=1)
                    for x in np.nonzero(ok)[0]:
                        red=len(AN)-len(set(R[x].tolist()))
                        if red>=6: out.append((pa_n,ca_n,s,kind,Per,red,int(x)))
    return out
if __name__=='__main__':
    # positive control
    random.seed(1);pt=[random.choice(ST) for _ in range(97)]
    for i,c in ANCH.items(): pt[i]=c
    key=[ST.index(c) for c in "MEDUSA"];M=colperms(97,7);m=M[1234];ct=['']*97
    for i in range(97): ct[m[i]]=KA[(ST.index(pt[i])+key[i%6])%26]
    print('control:',[h for h in scan(''.join(ct),M) if h[-1]==1234][:2])
    I=np.arange(97,dtype=np.int16)[None,:]
    print('identity:',scan(CT0,I))
    for W in range(2,9):
        h=scan(CT0,colperms(97,W));print('W',W,'hits',len(h),h[:5])
