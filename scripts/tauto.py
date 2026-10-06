import itertools,sys,random,numpy as np
from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
N=97;AN=np.array(sorted(ANCH));
def colmap(W,order):
    L=N//W;r=N%W;lens=[L+(1 if c<r else 0) for c in range(W)];idx=[]
    for c in order: idx+=[row*W+c for row in range(lens[c])]
    return idx            # Y[n]=X[idx[n]]
def perms(W):
    M=np.array([colmap(W,o) for o in itertools.permutations(range(W))],dtype=np.int16)
    Minv=np.argsort(M,axis=1).astype(np.int16)
    return M,Minv
THR_PT=8
def run(ct,Ws,label=''):
    hits=[]
    for W in Ws:
        M,Minv=perms(W)
        for dname,m,minv in (('enc',M,Minv),('dec',Minv,M)):
            # m: Y[n]=X[m[n]], minv: X position i -> Y position minv[i]
            for an,A in (('ST',ST),('KA',KA)):
                ci=np.array([A.index(x) for x in ct]);P=np.array([A.index(ANCH[i]) for i in AN])
                Pfull=np.full(N,-1);Pfull[AN]=P
                for rule in (0,1,2):
                    def key(c,p): return (c-p)%26 if rule==0 else (c+p)%26 if rule==1 else (p-c)%26
                    # ---- order A: autokey on PT order, then transposition. X[i]=ct[minv[i]]
                    X=ci[minv]                      # (np, 97)
                    kA=key(X[:,AN],P)               # (np,24)
                    for L in range(1,N):
                        sel=AN>=L
                        if sel.sum()>=20:
                            prev=X[:,AN[sel]-L]
                            for sg in (1,-1):
                                d=(kA[:,sel]-sg*prev)%26
                                ok=np.all(d==d[:,:1],axis=1)
                                for r in np.nonzero(ok)[0]: hits.append(('A-CTauto',W,dname,an,rule,L,sg,r))
                        # PT autokey order A: perm-independent pairs
                        pr=[(a,np.searchsorted(AN,AN[a]-L)) for a in range(len(AN)) if (AN[a]-L) in ANCH]
                        if len(pr)>=THR_PT:
                            a_=np.array([x for x,_ in pr]);b_=np.array([y for _,y in pr])
                            for sg in (1,-1):
                                d=(kA[:,a_]-sg*P[b_])%26
                                ok=np.all(d==d[:,:1],axis=1)
                                for r in np.nonzero(ok)[0]: hits.append(('A-PTauto',W,dname,an,rule,L,sg,r))
                    # ---- order B: transposition then autokey on Y order. anchor i sits at Y pos j=minv[i]
                    J=minv[:,AN].astype(int)            # (np,24)
                    cj=ci[J];kB=key(cj,P)
                    Yk=np.full(M.shape,-1,dtype=np.int16);np.put_along_axis(Yk,J,P[None,:].repeat(len(J),0),axis=1)
                    for L in range(1,N):
                        JL=J-L;valid=JL>=0
                        cprev=ci[np.clip(JL,0,None)]
                        for sg in (1,-1):
                            d=(kB-sg*cprev)%26
                            d=np.where(valid,d,-1)
                            # all valid equal, count>=20
                            ref=np.max(d,axis=1)
                            ok=(np.sum(valid,axis=1)>=20)&np.all((d==ref[:,None])|(~valid),axis=1)
                            for r in np.nonzero(ok)[0]: hits.append(('B-CTauto',W,dname,an,rule,L,sg,r))
                        # PT autokey order B: key_j = Y_{j-L}
                        yp=np.take_along_axis(Yk,np.clip(JL,0,None),axis=1).astype(int);v=valid&(yp>=0)
                        cnt=v.sum(axis=1);cand=np.nonzero(cnt>=THR_PT)[0]
                        if len(cand):
                            for sg in (1,-1):
                                d=(kB[cand]-sg*yp[cand])%26;vv=v[cand]
                                d=np.where(vv,d,-1);ref=d.max(axis=1)
                                ok=np.all((d==ref[:,None])|(~vv),axis=1)
                                for r in cand[np.nonzero(ok)[0]]: hits.append(('B-PTauto',W,dname,an,rule,L,sg,r))
        print(label,'W',W,'hits so far',len(hits),flush=True)
    return hits
if __name__=='__main__':
    mode=sys.argv[1]
    if mode=='ctrl':
        src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
        pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
        from autokey import enc_ctauto,enc_ptauto
        order=(3,0,5,1,4,2);idx=colmap(6,order)
        # A: autokey then transpose
        x=enc_ctauto(pt,KA,5,'CLOCK',0);ctA=''.join(x[i] for i in idx)
        h=run(ctA,[6],'ctrlA');print([t for t in h if t[0]=='A-CTauto'][:3])
        # B: transpose then PT-autokey
        y=''.join(pt[i] for i in idx);ctB=enc_ptauto(y,ST,3,'BER',1)
        h=run(ctB,[6],'ctrlB');print([t for t in h if t[0]=='B-PTauto'][:3])
    elif mode=='null':
        random.seed(int(sys.argv[2]));s=list(CT0);random.shuffle(s)
        h=run(''.join(s),range(2,9),'null');print('NULL hits',len(h),h[:5])
    else:
        h=run(CT0,[int(w) for w in sys.argv[2].split(',')],'K4');print('K4 hits',len(h),h[:20])
