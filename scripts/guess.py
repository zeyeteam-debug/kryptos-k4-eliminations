import numpy as np,random,itertools
from gen import run
from stack import CT0,ANCH
TWO='OF TO IS IN AT ON BY AS OR IT AN BE UP SO NO WE HE GO DO US IF AM XX'.split()
FOUR=('FROM NEAR PAST OVER INTO ONTO UPON WITH THIS THAT THEN FACE SEEK FIND VIEW LOOK SEES TOWN WEST EAST SIDE GATE WALL '
 'TIME HOUR WHEN WERE WILL BACK DOWN AWAY LEFT TURN READ TOLD SAID LIKE JUST ONLY ALSO BOTH EACH SOME MORE MOST MUCH EVEN '
 'HERE HAVE BEEN THEY THEM WHAT YOUR INTO ABLE AREA CITY DOOR ROOM MARK SIGN STEP WALK PATH ROAD LINE HALF NEXT LAST FIRST '
 'ATTH THEX XTHE ANDT OFTH TOTH ATTE WITH ALEX PLAT ZONE DIAL HAND FACE TOWR SHOW POIN').split()
FOUR=[w for w in FOUR if len(w)==4]+[a+b for a in TWO for b in TWO if a!='XX' and b!='XX']
FOUR=list(dict.fromkeys(FOUR))
Wp=[i for i,c in enumerate(CT0) if c=='W'];starts=[0]+[w+1 for w in Wp]
def segcheck(ct,anch):
    loc={i:i-max(s for s in starts if s<=i) for i in anch}
    s1=[i for i in anch if 21<=i<=35];s4=[i for i in anch if 59<=i<=73]
    pairs=[(a,b) for a in s1 for b in s4 if loc[a]==loc[b]]
    bad=sum((anch[a]==anch[b])!=(ct[a]==ct[b]) for a,b in pairs);strong=sum(anch[a]==anch[b] for a,b in pairs)
    return bad,strong,len(pairs)
def scan(ct):
    res=[]
    one=np.arange(97,dtype=np.int16)[None,:]
    for t in TWO:
        for f in FOUR:
            a=dict(ANCH);a[34]=t[0];a[35]=t[1]
            for j,ch in enumerate(f): a[59+j]=ch
            h=run(ct,a,one,[(t,f)],Pmax=24,RTHR=4,ANYSTRONG=3)
            bad,strong,n=segcheck(ct,a)
            res.append((t,f,len(h),h[:2],bad,strong,n))
    return res
if __name__=='__main__':
    print('combinations',len(TWO)*len(FOUR))
    r=scan(CT0)
    print('K4 periodic-key hits:',[(x[0],x[1],x[3]) for x in r if x[2]][:10])
    seg=[(x[0],x[1],x[5]) for x in r if x[4]==0];print('K4 W-segment (any alphabet) survivors:',len(seg),'with strong>=2:',[s for s in seg if s[2]>=2][:20])
    # null
    cnt=[];cnt2=[]
    for sd in range(10):
        l=list(CT0);rng=random.Random(sd)
        nonW=[i for i,c in enumerate(CT0) if c!='W'];lt=[CT0[i] for i in nonW];rng.shuffle(lt)
        for i,c in zip(nonW,lt): l[i]=c
        rr=scan(''.join(l));cnt.append(sum(1 for x in rr if x[2]));cnt2.append((sum(1 for x in rr if x[4]==0),sum(1 for x in rr if x[4]==0 and x[5]>=2)))
    print('null periodic hits per shuffle:',cnt);print('null W-segment survivors (all, strong>=2):',cnt2)
