# Paper fortune teller (cootie catcher) fold transposition + flower/lotus/fortune-teller vocabulary alphabets
import numpy as np,itertools,random,time
from stack import CT0,ANCH,KA,ST
from gen import run
from qii_more import scan as qii
from flower import vogel_perms
N=97;n=10;c=n/2
def fold_states(flip=True,half=True,blintz2=True):
    X=np.array([j+.5 for i in range(n) for j in range(n)]);Y=np.array([i+.5 for i in range(n) for j in range(n)]);Z=np.zeros(n*n)
    def fold(mask,fx,fy):
        nonlocal X,Y,Z
        zm=Z.max();x,y=X[mask].copy(),Y[mask].copy();X[mask],Y[mask]=fx(x,y),fy(x,y);Z[mask]=2*zm+1-Z[mask]
    # blintz 1: four corners to centre
    fold(X+Y<c,lambda x,y:c-y,lambda x,y:c-x)
    fold((n-X)+Y<c,lambda x,y:y+c,lambda x,y:x-c)
    fold(X+(n-Y)<c,lambda x,y:y-c,lambda x,y:x+c)
    fold(X+Y>n+c,lambda x,y:n+c-y,lambda x,y:n+c-x)
    if flip: X=n-X;Z=-Z
    if blintz2:
        q=c/2
        fold(Y<q,lambda x,y:x,lambda x,y:2*q-y)
        fold(Y>n-q,lambda x,y:x,lambda x,y:2*(n-q)-y)
        fold(X<q,lambda x,y:2*q-x,lambda x,y:y)
        fold(X>n-q,lambda x,y:2*(n-q)-x,lambda x,y:y)
    if half: fold(Y>c,lambda x,y:x,lambda x,y:n-y)
    return X,Y,Z
def writes():
    W={}
    W['rows']=list(range(n*n));W['cols']=[i*n+j for j in range(n) for i in range(n)]
    W['snake']=[i*n+(j if i%2==0 else n-1-j) for i in range(n) for j in range(n)]
    W['spiral']=None
    # clockwise spiral from top-left
    g=[[i*n+j for j in range(n)] for i in range(n)];s=[]
    while g: s+=g.pop(0);g=[list(r) for r in zip(*g)][::-1]
    W['spiral']=s
    return W
def cootie_perms():
    P=[];L=[]
    for flip,half,b2 in itertools.product((True,False),repeat=3):
        X,Y,Z=fold_states(flip,half,b2);Xr,Yr,Zr=np.round(X,3),np.round(Y,3),Z
        reads={'pos_rows_zup':lambda k:(Yr[k],Xr[k],Zr[k]),'pos_rows_zdn':lambda k:(Yr[k],Xr[k],-Zr[k]),
               'pos_cols_zup':lambda k:(Xr[k],Yr[k],Zr[k]),'pos_cols_zdn':lambda k:(Xr[k],Yr[k],-Zr[k]),
               'z_up_rows':lambda k:(Zr[k],Yr[k],Xr[k]),'z_dn_rows':lambda k:(-Zr[k],Yr[k],Xr[k]),
               'z_up_cols':lambda k:(Zr[k],Xr[k],Yr[k]),'z_dn_cols':lambda k:(-Zr[k],Xr[k],Yr[k]),
               'angle_z':lambda k:(round(np.arctan2(Yr[k]-c,Xr[k]-c),4),Zr[k],k),'angle_zdn':lambda k:(round(np.arctan2(Yr[k]-c,Xr[k]-c),4),-Zr[k],k)}
        for wn,w in writes().items():
            for pad in ('end','start','corners'):
                padcells={'end':w[97:],'start':w[:3],'corners':None}[pad]
                if pad=='corners': cells=[x for x in w if x not in (0,n-1,n*n-1)][:97]
                else: cells=[x for x in w if x not in padcells]
                for rn,f in reads.items():
                    order=sorted(cells,key=lambda k:(f(k),k));rank={k:r for r,k in enumerate(order)}
                    m=np.array([rank[cells[a]] for a in range(N)])
                    P+=[m,np.argsort(m)];L+=[(flip,half,b2,wn,pad,rn,'enc'),(flip,half,b2,wn,pad,rn,'dec')]
    P=np.array(P,dtype=np.int16);u,ix=np.unique(P,axis=0,return_index=True);ix=sorted(ix)
    return P[ix],[L[i] for i in ix]
def kalpha(w):
    s=''
    for ch in w+ST:
        if ch not in s: s+=ch
    return s
WORDS=['LOTUS','BLUELOTUS','LOTOS','NELUMBO','NYMPHAEA','NYMPHAEACAERULEA','SESHEN','SESHESHEN','PAPYRUS','SUNFLOWER','HELIANTHUS',
 'ROSE','LILY','TULIP','IRIS','DAISY','POPPY','ORCHID','LAVENDER','HOROLOGIUMFLORAE','LINNAEUS','FLOWERCLOCK','BLUMENUHR','FLORALCLOCK',
 'PHYLLOTAXIS','FIBONACCI','GOLDENANGLE','PETAL','PETALS','BLOOM','BLOSSOM','MORNINGGLORY','FOUROCLOCK','MIRABILIS','DANDELION','CHICORY',
 'NIGHTBLOOMINGCEREUS','SHOSHAN','ANTHOS','LOULOUDI','ROSETTE','ORIGAMI','COOTIECATCHER','FORTUNETELLER','PAPERFORTUNETELLER','CHATTERBOX',
 'SALTCELLAR','WHIRLYBIRD','FLOWER','OPENCLOSE','KRYPTOSLOTUS','LOTUSKRYPTOS']
def alph_scan(ct,M,pairs,Pmax=26):
    AN=np.array(sorted(ANCH));PR=np.array(list(itertools.combinations(range(len(AN)),2)));hits=[]
    for lab,pa,ca in pairs:
        C=np.array([ca.index(x) for x in ct])[M[:,AN]];Pp=np.array([pa.index(ANCH[i]) for i in AN])
        for s in (1,-1):
            K=(C-s*Pp)%26;same=K[:,PR[:,0]]==K[:,PR[:,1]]
            for Per in range(1,Pmax+1):
                for rn,R in (('pt',np.broadcast_to(AN%Per,M[:,AN].shape)),('ct',M[:,AN]%Per)):
                    need=R[:,PR[:,0]]==R[:,PR[:,1]];ok=np.all(same|~need,axis=1)
                    for x in np.nonzero(ok)[0]:
                        if len(AN)-len(set(R[x].tolist()))>=6: hits.append((lab,s,Per,rn,int(x)))
    return hits
if __name__=='__main__':
    t=time.time();P,L=cootie_perms();print('cootie fold perms',len(P),flush=True)
    SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for a,ch in ANCH.items(): pt[a]=ch
    i=len(P)//3;m=P[i];ct=['?']*N
    for a in range(N): ct[m[a]]=KA[(KA.index(pt[a])+KA.index('FLOWERS'[a%7]))%26]
    print('fold control',[h for h in run(''.join(ct),ANCH,P[i:i+1],[L[i]],Pmax=12) if h[0]=='key'][:1],flush=True)
    for nm,x in (('K4',CT0),('null',SH)):
        H=run(x,ANCH,P,L,Pmax=24);q=qii(x,P,24)
        print(nm,'fold: periodic hits',len(H),H[:3],'| QI/II mixed',q,flush=True)
    # flower-keyword alphabets, Quagmire I-IV, with identity + fold + flower perms
    FP,FL=vogel_perms();ALLP=np.vstack([np.arange(N,dtype=np.int16)[None],P,FP])
    pairs=[]
    for w in WORDS:
        A=kalpha(w)
        for lab,pa,ca in ((w+':QIII',A,A),(w+':QI',A,ST),(w+':QII',ST,A),(w+':W/KA',A,KA),(w+':KA/W',KA,A)):
            pairs.append((lab,pa,ca))
    # control: plant QIII with LOTUS alphabet, period 5, identity
    A=kalpha('LOTUS');ct=''.join(A[(A.index(pt[a])+A.index('BLOOM'[a%5]))%26] for a in range(N))
    print('alpha control',[h for h in alph_scan(ct,ALLP[:1],[('LOTUS:QIII',A,A)]) if h[2]==5][:1],flush=True)
    for nm,x in (('K4',CT0),('null',SH)):
        H=[]
        for c0 in range(0,len(ALLP),4000): H+=[(h[0],h[1],h[2],h[3],h[4]+c0) for h in alph_scan(x,ALLP[c0:c0+4000],pairs)]
        print(nm,'flower alphabets x',len(pairs),'pairs x',len(ALLP),'perms: hits',len(H),H[:5],round(time.time()-t),flush=True)
