import numpy as np,random,math
exec(open('board.py').read().split('# K4 cell')[0])
from stack import CT0,ANCH
from sculpt import TAB
KA=KRY;ST=STD;AN=np.array(sorted(ANCH))
def gridarr(G):
    R=len(G);W=max(len(r) for r in G);a=np.full((R,W),-1)
    return G,R,W
def needs(ct):
    out={}
    for an,A in (('KA',KA),('ST',ST)):
        c=np.array([A.index(ct[i]) for i in AN]);p=np.array([A.index(ANCH[i]) for i in AN])
        for r in range(3): out[(an,r)]=(A,[(c-p)%26,(c+p)%26,(p-c)%26][r])
    return out
def scan(ct,want=None):
    hits=[];nd=needs(ct);best=0
    for gname,G in (('tableau',TAB),('panel',LINES),('tableau_rev',TAB[::-1]),('panel_rev',LINES[::-1])):
        R=len(G);Wc=min(len(r) for r in G)
        lut={key:np.array([[ (A.index(ch) if ch in A else -1) for ch in row.ljust(40)] for row in G]) for key,(A,_) in nd.items()}
        for Amp in range(1,7):
            for lam in (8,10,12,14,16,20,24,28,32,40,48,64,97):
                for phi in range(lam):
                    col=AN%Wc;pas=AN//Wc
                    wav=np.rint(Amp*np.sin(2*math.pi*(col+phi)/lam)).astype(int)
                    for dlt in range(1,2*Amp+2):
                        for r0 in range(Amp,R):
                            row=r0+pas*dlt+wav
                            if row.max()>=R or row.min()<0: continue
                            for key,(A,need) in nd.items():
                                kv=lut[key][row,col]
                                m=int((kv==need).sum());best=max(best,m)
                                if m==24: hits.append((gname,Amp,lam,phi,dlt,r0)+key)
    return hits,best
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    Amp,lam,phi,dlt,r0=3,16,5,4,4;Wc=min(len(r) for r in TAB);key=[]
    for i in range(97):
        col=i%Wc;row=r0+(i//Wc)*dlt+round(Amp*math.sin(2*math.pi*(col+phi)/lam));key.append(TAB[row][col])
    ctrl=''.join(KA[(KA.index(p)+KA.index(k))%26] for p,k in zip(pt,key))
    h,b=scan(ctrl);print('control:',[x for x in h if x[:6]==('tableau',3,16,5,4,4)][:1],'total',len(h),flush=True)
    h,b=scan(CT0);print('K4 full:',h[:5],'best',b,'/24',flush=True)
    l=list(CT0);random.Random(4).shuffle(l);print('null best',scan(''.join(l))[1])
