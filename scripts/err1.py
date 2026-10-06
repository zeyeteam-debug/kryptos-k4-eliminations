import numpy as np,random,sys
from gen import run,colperms
from stack import CT0,ANCH
def variants(ct,anch):
    yield ('exact',),ct,anch
    for i in sorted(anch):
        a=dict(anch);del a[i];yield ('drop',i),ct,a
    for j in range(len(ct)):
        yield ('del',j),ct[:j]+ct[j+1:],anch
    for j in range(len(ct)+1):
        yield ('ins',j),ct[:j]+'?'+ct[j:],anch
def direct(ct,anch):
    h=[]
    for lab,c,a in variants(ct,anch):
        M=np.arange(len(c),dtype=np.int16)[None,:]
        h+=[x[:-1]+(lab,) for x in run(c,a,M,[lab],Pmax=24,RTHR=4,ANYSTRONG=3)]
    return h
def columnar(ct,anch,Ws,indel=True):
    h=[];cache={}
    for lab,c,a in variants(ct,anch):
        if not indel and lab[0] in ('del','ins'): continue
        for W in Ws:
            key=(len(c),W)
            if key not in cache: cache[key]=colperms(len(c),W)
            Ms=cache[key]
            h+=[x[:-1]+(lab,W,x[-1]) for x in run(c,a,Ms,list(range(len(Ms))),Pmax=12,RTHR=6)]
    return h
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
    # control 1: periodic KA key P=7, then one ciphertext letter corrupted at anchor 66 and one extra letter inserted at 50
    ct=''.join(KA[(KA.index(p)+KA.index('KRYPTOS'[i%7]))%26] for i,p in enumerate(pt))
    c1=list(ct);c1[66]='Z' if c1[66]!='Z' else 'Y';c1=''.join(c1)
    c2=ct[:50]+'Q'+ct[50:]
    for nm,c in (('typo',c1),('extra letter',c2)):
        h=direct(c,ANCH);print('control',nm,':',[x for x in h if x[0]=='key' and x[1]=='KA' and x[4]==7][:2],'all hits',len(h))
    nh=[]
    for sd in range(10):
        l=list(CT0);random.Random(sd).shuffle(l);nh.append(len(direct(''.join(l),ANCH)))
    print('null direct hits per shuffle:',nh)
    h=direct(CT0,ANCH);print('K4 direct hits',len(h),h[:15],flush=True)
    if len(sys.argv)>1:
        Ws=[int(x) for x in sys.argv[1].split(',')]
        from tauto import colmap
        idx=colmap(6,(2,0,5,3,1,4));x=''.join(KA[(KA.index(p)+KA.index('CLOCK'[i%5]))%26] for i,p in enumerate(pt))
        cc=''.join(x[i] for i in idx);cc=cc[:30]+'R'+cc[31:] if cc[30]!='R' else cc[:30]+'S'+cc[31:]
        hc=columnar(cc,ANCH,[6],indel=False);print('control columnar W6 + typo:',len(hc),hc[:2],flush=True)
        l=list(CT0);random.Random(99).shuffle(l);print('null columnar hits',len(columnar(''.join(l),ANCH,Ws,indel=False)),flush=True)
        hk=columnar(CT0,ANCH,Ws,indel=False);print('K4 columnar(drop-one) hits',len(hk),hk[:15])
