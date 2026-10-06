import numpy as np,random,sys
from stackT import check,ST,KA
from stack import CT0,ANCH
from gen import colperms
AN=np.array(sorted(ANCH))
def scan(ct,Ms,labs,amax=12,thr=4):
    hits=[]
    for an,A in (('ST',ST),('KA',KA)):
        C=np.array([A.index(x) for x in ct]);P=np.array([A.index(ANCH[i]) for i in AN])
        for x in range(len(Ms)):
            q=Ms[x][AN]
            for e in (1,-1):
                D=((C[q]-e*P)%26).tolist()
                for a in range(2,amax+1):
                    I=(q%a).tolist()
                    for b in range(a+1,amax+1):
                        if b%a==0: continue
                        r=check(I,(q%b).tolist(),D,a,b)
                        if r>=thr: hits.append((an,e,a,b,r,labs[x]))
    return hits
if __name__=='__main__':
    from tauto import colmap
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    idx=colmap(5,(3,1,4,0,2));y=''.join(pt[i] for i in idx)
    ctrl=''.join(KA[(KA.index(c)+KA.index('OCHRE'[j%5])+KA.index('VERDIGRIS'[j%9]))%26] for j,c in enumerate(y))
    M5=colperms(97,5);print('control T5->Q(OCHRE)->Q(VERDIGRIS):',[h for h in scan(ctrl,M5,list(range(len(M5)))) if h[2:4]==(5,9)][:2],flush=True)
    for W in [int(w) for w in sys.argv[1].split(',')]:
        Ms=colperms(97,W);labs=[(W,i) for i in range(len(Ms))]
        l=list(CT0);random.Random(W).shuffle(l)
        print('W',W,'null',len(scan(''.join(l),Ms,labs)),'K4',scan(CT0,Ms,labs)[:5],flush=True)
