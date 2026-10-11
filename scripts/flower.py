# Flower hypotheses: phyllotaxis (Vogel sunflower spiral) reading orders, Fibonacci decimations, opening/closing readings;
# flower-clock / lotus vocabulary as keys.
import numpy as np,math,random,time
from stack import CT0,ANCH
from gen import run
from qii_more import scan as qii
N=97
def vogel_perms():
    P=[];L=[]
    for name,ang in (('golden',137.50776),('golden2',222.49224),('lucas',99.50192),('noble3',151.1374)):
        th=[(n*ang)%360 for n in range(N)];r=[math.sqrt(n) for n in range(N)]
        orders={}
        orders['by_angle']=sorted(range(N),key=lambda n:(th[n],r[n]))
        orders['by_angle_desc']=sorted(range(N),key=lambda n:(-th[n],r[n]))
        orders['by_angle_out_in']=sorted(range(N),key=lambda n:(th[n],-r[n]))
        for F in (3,5,8,13,21,34,55,89):          # parastichy families: n, n+F, n+2F,... read family by family
            orders[f'parastichy{F}']=[n for s in range(F) for n in range(s,N,F)]
            orders[f'parastichy{F}_alt']=[n for s in range(F) for n in (list(range(s,N,F)) if s%2==0 else list(range(s,N,F))[::-1])]  # open/close
        for k in (4,6,8,12,24):                   # radial rings (petal layers), opening = outward, closing = inward
            ring=[int(r[n]*k/ math.sqrt(N)) for n in range(N)]
            orders[f'rings{k}_out']=sorted(range(N),key=lambda n:(ring[n],th[n]))
            orders[f'rings{k}_alt']=sorted(range(N),key=lambda n:(ring[n],th[n] if ring[n]%2==0 else -th[n]))
        for on,o in orders.items():
            o=np.array(o);P+=[o,np.argsort(o)];L+=[(name,on,'a'),(name,on,'b')]
    for k in range(2,N):                          # decimation by k (Fibonacci = 13,21,34,55,89) any start
        for st in range(N):
            o=np.array([(st+k*j)%N for j in range(N)]);P+=[o];L+=[('decim',k,st)]
    P=np.array(P,dtype=np.int16);u,ix=np.unique(P,axis=0,return_index=True);ix=sorted(ix)
    return P[ix],[L[i] for i in ix]
if __name__=='__main__':
    t=time.time();P,L=vogel_perms();print('flower/decimation perms',len(P),flush=True)
    SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    # control: golden parastichy13 + KA vig period 6
    KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';i=[x for x,l in enumerate(L) if l[:2]==('golden','parastichy13')][0];m=P[i]
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for a,c in ANCH.items(): pt[a]=c
    ct=['?']*N
    for a in range(N): ct[m[a]]=KA[(KA.index(pt[a])+KA.index('LOTUSX'[a%6]))%26]
    print('control',[h for h in run(''.join(ct),ANCH,P[i:i+1],[L[i]],Pmax=12) if h[0]=='key'][:1],flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        H=[];q=0
        for c in range(0,len(P),20000):
            H+=run(ct,ANCH,P[c:c+20000],L[c:c+20000],Pmax=24);q+=qii(ct,P[c:c+20000],12)
        print(nm,'periodic ST/KA/any hits',len(H),[h for h in H if h[0]=='key'][:3],'any',sum(1 for h in H if h[0]=='any'),'| Quagmire I/II mixed',q,round(time.time()-t),flush=True)
