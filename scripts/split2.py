# K4 = two segments with independent systems. Split s in 34..63: segment A=[0,s) holds EASTNORTHEAST, B=[s,97) holds BERLINCLOCK.
# Each segment: identity or its own columnar transposition (widths 2..7, all orders, both directions) + periodic key
# (ST/KA Vig/Beau, periods 1..12, pt- or ct-indexed; and arbitrary alphabets). Count consistent configs per segment vs nulls.
import numpy as np,random,time,sys
from stack import CT0,ANCH
from gen import run,colperms
def seg_perms(L):
    P=[np.arange(L,dtype=np.int16)];lab=[('id',)]
    for W in range(2,8):
        M=colperms(L,W);P+=list(M);lab+=[('col',W,k) for k in range(len(M))]
    return np.array(P,dtype=np.int16),lab
def run_split(ct):
    res={}
    cache={}
    for s in range(34,64):
        for name,(a,b) in (('A',(0,s)),('B',(s,97))):
            L=b-a;sub=ct[a:b];an={i-a:c for i,c in ANCH.items() if a<=i<b}
            if L not in cache: cache[L]=seg_perms(L)
            P,lab=cache[L]
            H=run(sub,an,P,lab,Pmax=12,RTHR=6,ANYSTRONG=4)
            res[(s,name)]=len(H)
    return res
if __name__=='__main__':
    t=time.time();which=sys.argv[1]
    if which=='K4': ct=CT0
    elif which=='ctrl':
        # plant: split at 50; A: columnar 5 + KA vig period 4 (pt-indexed); B: identity + ST beau period 3
        KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
        pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
        for i,c in ANCH.items(): pt[i]=c
        PA,_=seg_perms(50);m=PA[200];ctA=['?']*50
        for i in range(50): ctA[m[i]]=KA[(KA.index(pt[i])+KA.index('LUCI'[i%4]))%26]
        ctB=''.join(ST[(ST.index('DOG'[(i-50)%3])-ST.index(pt[i]))%26] for i in range(50,97))
        ct=''.join(ctA)+ctB
    else:
        l=list(CT0);random.Random(int(which[4:])).shuffle(l);ct=''.join(l)
    r=run_split(ct)
    both=[s for s in range(34,64) if r[(s,'A')] and r[(s,'B')]]
    print(which,'A-consistent splits',sum(1 for s in range(34,64) if r[(s,'A')]),'B-consistent',sum(1 for s in range(34,64) if r[(s,'B')]),
          'both',both,'| total configs A',sum(r[(s,'A')] for s in range(34,64)),'B',sum(r[(s,'B')] for s in range(34,64)),'| s=50:',r[(50,'A')],r[(50,'B')],round(time.time()-t),flush=True)
