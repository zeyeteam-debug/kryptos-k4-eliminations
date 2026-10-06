import numpy as np,random,sys
from gen import run,colperms
from stack import CT0,ANCH
def reduce(ct,nullset):
    keep=[i for i,c in enumerate(ct) if c not in nullset]
    newpos={old:n for n,old in enumerate(keep)}
    anch={newpos[i]:p for i,p in ANCH.items() if i in newpos}
    return ''.join(ct[i] for i in keep),anch
def test(ct,anch,Ws):
    h=[];L=len(ct)
    h+=[x for x in run(ct,anch,np.arange(L,dtype=np.int16)[None,:],[('direct',)],Pmax=24,RTHR=4,ANYSTRONG=3)]
    for W in Ws:
        Ms=colperms(L,W);h+=run(ct,anch,Ms,[('col',W,i) for i in range(len(Ms))])
    return h
if __name__=='__main__':
    KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
    # control: encrypt with KA key period 6, then insert W nulls at K4's W positions
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:92]);pt[20:33]='EASTNORTHEAST';pt[59:70]='BERLINCLOCK';pt=''.join(pt)
    c=[KA[(KA.index(p)+KA.index('SHADOW'[i%6]))%26] for i,p in enumerate(pt)]
    c=[x if x!='W' else 'V' for x in c]          # keep W only as null
    for pos in [20,36,48,58,74]: c.insert(pos,'W')
    c=''.join(c);r,a=reduce(c,'W');print('control found:',[x for x in test(r,a,[]) if x[0]=='key' and x[4]==6][:1])
    Ws=[int(x) for x in sys.argv[1].split(',')] if len(sys.argv)>1 else []
    r,a=reduce(CT0,'W');print('W as nulls: text',len(r),'anchors',len(a));h=test(r,a,Ws);print('K4 W-null hits',len(h),h[:10],flush=True)
    # every single letter as null, direct only, for comparison
    for L in 'ABCDEFGHIJKLMNOPQRSTUVXYZ':
        r,a=reduce(CT0,L);h=test(r,a,[])
        if h: print('letter',L,'as null ->',h[:3])
    print('done single letters')
