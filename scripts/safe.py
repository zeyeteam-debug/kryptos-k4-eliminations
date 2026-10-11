# Safe wheel-pack keystream: key is piecewise linear in ct position, 3 phases with slopes s1,s2,s3 in {+-1,+-2,+-3}
# (mechanism: s=(-3,+2,-1) or (+3,-2,+1) = all wheels / two wheels / one wheel turning), any phase boundaries,
# on top of transpositions (identity, columnar 2-8, alternating dial, K3-style double rotation).
import numpy as np,itertools,random,time
from stack import CT0,ANCH
from gen import colperms
from mech import dial_perms,dial
from k3rot import compose
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';AN=np.array(sorted(ANCH));N=97
def perms():
    P=[np.arange(N)];L=[('id',)]
    for W in range(2,9):
        M=colperms(N,W);P+=list(M);L+=[('col',W,k) for k in range(len(M))]
    D,DL=dial_perms([(1,1,97),(2,1,97),(3,1,26),(4,1,12)]);P+=list(D);L+=DL
    for w1 in range(2,97):
        for w2 in range(2,97):
            for d1 in (True,False):
                for d2 in (True,False):
                    idx=compose(N,w1,d1,w2,d2);inv=np.empty(N,int);inv[idx]=np.arange(N)
                    P+=[inv,idx];L+=[('rot',w1,d1,w2,d2,'e'),('rot',w1,d1,w2,d2,'d')]
    return np.array(P,dtype=np.int16),L
SL=[-3,-2,-1,1,2,3]
SEQS=list(itertools.product(SL,repeat=3))
MECH={(-3,2,-1),(3,-2,1)}
def search(ct,Ms,labs,chunk=40000):
    hits=[]
    for an,A in (('ST',ST),('KA',KA)):
        lut=np.zeros(128,int);lut[[ord(c) for c in A]]=np.arange(26);cl=np.array([ord(c) for c in ct]);P=np.array([A.index(ANCH[i]) for i in AN])
        for rule in (0,1):
            for c0 in range(0,len(Ms),chunk):
                M=Ms[c0:c0+chunk].astype(np.int64);pos=M[:,AN];C=lut[cl[pos]]
                K=(C-P)%26 if rule==0 else (C+P)%26
                o=np.argsort(pos,axis=1);ps=np.take_along_axis(pos,o,1);Ks=np.take_along_axis(K,o,1)
                EQ={s:((Ks-s*ps)%26)[:,1:]==((Ks-s*ps)%26)[:,:-1] for s in SL}   # n x 23
                PRE={s:np.argmin(np.c_[EQ[s],np.zeros(len(M),bool)],axis=1) for s in SL}   # leading trues
                SUF={s:np.argmin(np.c_[EQ[s][:,::-1],np.zeros(len(M),bool)],axis=1) for s in SL}
                for s1,s2,s3 in SEQS:
                    x=PRE[s1]+1; y=23-SUF[s3]          # segment1 = crib idx [0,x), segment3 = [y,24)
                    ok=np.zeros(len(M),bool)
                    full=(x>=y)                          # segments 1 and 3 already cover / overlap
                    ok|=full
                    rest=np.nonzero(~full)[0]
                    if len(rest):
                        e=EQ[s2][rest];xr=x[rest];yr=y[rest]
                        cs=np.c_[np.zeros(len(rest),int),np.cumsum(~e,axis=1)]
                        # middle segment crib idx [xr,yr): needs eq[xr .. yr-2] all true
                        lo=xr;hi=np.maximum(yr-1,xr)
                        bad=cs[np.arange(len(rest)),hi]-cs[np.arange(len(rest)),lo]
                        ok[rest]=bad==0
                    for j in np.nonzero(ok)[0]: hits.append(((s1,s2,s3) in MECH,(s1,s2,s3),an,rule,labs[c0+j]))
    return hits
if __name__=='__main__':
    t=time.time();Ms,labs=perms();print('perms',len(Ms),round(time.time()-t),flush=True)
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    d=np.array(dial((15,25,2),'R'));key=[(5-3*p if p<30 else (5-90+2*(p-30) if p<70 else 9-1*(p-70)))%26 for p in range(97)]
    ct=['?']*97
    for i,c in enumerate(pt): ct[d[i]]=KA[(KA.index(c)+key[d[i]])%26]
    h=search(''.join(ct),d[None,:].astype(np.int16),['ctrl']);print('control',[x for x in h if x[1]==(-3,2,-1)][:1],flush=True)
    SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    for nm,ct in (('K4',CT0),('null',SH)):
        h=search(ct,Ms,labs);print(nm,'hits',len(h),'mechanism-slope hits',sum(1 for x in h if x[0]),h[:5],round(time.time()-t),flush=True)
