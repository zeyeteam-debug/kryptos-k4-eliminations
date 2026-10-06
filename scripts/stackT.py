"""Sub(period a, PT index) -> Transposition -> Sub(period b, CT index), additive in one alphabet.
c[m[i]] - e*p[i] = K1[i%a] + K2[m[i]%b]  (e=+1/-1 covers Vig/Beaufort combos). Exact union-find with potentials mod 26."""
import numpy as np,sys,random,itertools
from stack import CT0,ANCH
from gen import colperms
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
AN=sorted(ANCH)
def check(I,J,D,na,nb):
    # nodes 0..na-1 = x, na..na+nb-1 = z ; constraint pot[x]-pot[z]=d
    par=list(range(na+nb));pot=[0]*(na+nb)
    def find(u):
        if par[u]==u: return u,0
        r,p=find(par[u]);par[u]=r;pot[u]=(pot[u]+p)%26;return r,pot[u]
    used=set();comps=0
    for u,v,d in zip(I,J,D):
        v=na+v;used.add(u);used.add(v)
        ru,pu=find(u);rv,pv=find(v)
        if ru==rv:
            if (pu-pv)%26!=d: return -1
        else:
            par[ru]=rv;pot[ru]=(pv+d-pu)%26
    roots=set(find(u)[0] for u in used)
    return len(I)-(len(used)-len(roots))
def scan(ct,Ms,labs,amax=8,bmax=8,thr=4):
    hits=[]
    for an,A in (('ST',ST),('KA',KA)):
        C=np.array([A.index(x) for x in ct]);P=np.array([A.index(ANCH[i]) for i in AN]);ANa=np.array(AN)
        for x in range(len(Ms)):
            q=Ms[x][ANa]
            for e in (1,-1):
                D=((C[q]-e*P)%26).tolist()
                for a in range(1,amax+1):
                    I=(ANa%a).tolist()
                    for b in range(2,bmax+1):
                        if a==1: continue
                        J=(q%b).tolist()
                        r=check(I,J,D,a,b)
                        if r>=thr: hits.append((an,e,a,b,r,labs[x]))
    return hits
if __name__=='__main__':
    Ws=[int(w) for w in sys.argv[1].split(',')]
    from tauto import colmap
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    k1='ABSCISSA';k2='PALIMPS'
    x=''.join(KA[(KA.index(p)+KA.index(k1[i%8]))%26] for i,p in enumerate(pt))
    idx=colmap(5,(3,1,4,0,2));y=''.join(x[i] for i in idx)
    ctrl=''.join(KA[(KA.index(c)+KA.index(k2[j%7]))%26] for j,c in enumerate(y))
    M5=colperms(97,5);h=scan(ctrl,M5,list(range(len(M5))));print('control found (a=8,b=7):',[t for t in h if t[2]==8 and t[3]==7][:2],'all',len(h),flush=True)
    for W in Ws:
        Ms=colperms(97,W);labs=[(W,i) for i in range(len(Ms))]
        l=list(CT0);random.Random(W).shuffle(l)
        print('W',W,'null hits',len(scan(''.join(l),Ms,labs)),'| K4 hits',scan(CT0,Ms,labs)[:10],flush=True)
