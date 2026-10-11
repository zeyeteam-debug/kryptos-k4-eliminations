# Quagmire I/II with fixed mixed pairs (pt ST / ct KA and pt KA / ct ST) over all our transposition families
import numpy as np,itertools,random,time
from stack import CT0,ANCH,KA,ST
from k3rot import compose
from mech import dial_perms
import dictsearch as D
AN=np.array(sorted(ANCH));PR=np.array(list(itertools.combinations(range(len(AN)),2)))
def scan(ct,M,Pmax=24):
    out=0
    for pa,ca in ((ST,KA),(KA,ST)):
        C=np.array([ca.index(c) for c in ct])[M[:,AN]];P=np.array([pa.index(ANCH[i]) for i in AN])
        for s in (1,-1):
            K=(C-s*P)%26;same=K[:,PR[:,0]]==K[:,PR[:,1]]
            for Per in range(1,Pmax+1):
                for R in (np.broadcast_to(AN%Per,M[:,AN].shape),M[:,AN]%Per):
                    need=R[:,PR[:,0]]==R[:,PR[:,1]];ok=np.all(same|~need,axis=1)
                    for x in np.nonzero(ok)[0]:
                        if len(AN)-len(set(R[x].tolist()))>=6: out+=1
    return out
def fams():
    N=97;R=[]
    for w1 in range(2,97):
        for w2 in range(2,97):
            for d1 in (True,False):
                for d2 in (True,False):
                    idx=compose(N,w1,d1,w2,d2);inv=np.empty(N,int);inv[idx]=np.arange(N);R+=[inv,idx]
    yield 'K3 rotations',np.array(R,dtype=np.int16)
    Pd,_=dial_perms([(1,1,97),(2,1,97),(3,1,26),(4,1,12)]);yield 'dial',Pd
    sp=[]
    for Wd in range(2,49):
        for corner in ('TL','TR','BL','BR'):
            for cw in (0,1):
                m=D.perm_from_read(D.spiral_order(97,Wd,corner,cw),97);sp+=[m,np.argsort(m)]
    yield 'spirals',np.array(sp,dtype=np.int16)
if __name__=='__main__':
    t=time.time();SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    I=np.arange(97,dtype=np.int16)[None,:];print('direct P<=24: K4',scan(CT0,I),'null',scan(SH,I),flush=True)
    for nm,M in fams():
        k=n=0
        for c in range(0,len(M),20000): k+=scan(CT0,M[c:c+20000],12);n+=scan(SH,M[c:c+20000],12)
        print(nm,len(M),'perms: K4',k,'null',n,round(time.time()-t),flush=True)
