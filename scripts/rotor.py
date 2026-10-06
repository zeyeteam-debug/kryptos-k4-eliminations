import random
from stack import ANCH,CT0
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
idx=sorted(ANCH)
def ok_rotor(ct,AL,N,dirn,sgn):
    # single rotor, unknown wiring R, rotor position s_i=sgn*floor(i/N); c = R(p+s)-s  (dirn 0) or p = R(c+s)-s (dirn 1)
    best=0
    for o in range(26):
        f={};g={};good=True;cons=0
        for i in idx:
            s=(sgn*(i//N)+o)%26;p=AL.index(ANCH[i]);c=AL.index(ct[i])
            x,y=((p+s)%26,(c+s)%26) if dirn==0 else ((c+s)%26,(p+s)%26)
            if x in f:
                if f[x]!=y: good=False;break
                cons+=1
            elif y in g: good=False;break
            else: f[x]=y;g[y]=x
        if good: return True
    return False
def run(ct):
    h=[]
    for an,AL in (('ST',ST),('KA',KA)):
        for N in range(1,49):
            for dirn in (0,1):
                for sgn in (1,-1):
                    if ok_rotor(ct,AL,N,dirn,sgn): h.append((an,N,dirn,sgn))
    return h
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    R=list(range(26));random.seed(5);random.shuffle(R)
    ct=''.join(ST[(R[(ST.index(ch)+i//3+7)%26]-(i//3+7))%26] for i,ch in enumerate(pt))
    print('control (N=3):',(('ST',3,0,1) in run(ct)))
    random.seed(9);fp=0
    for _ in range(100):
        s=list(CT0);random.shuffle(s);fp+=bool(run(''.join(s)))
    print('null shuffles with any hit: %d/100'%fp)
    print('K4:',run(CT0))
