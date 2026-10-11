# Literal combination-lock simulation. Dial ring = ST or KA letters. 4 phases:
#  P1 dir d,  wheels {1,2,3} move, T1 full turns then stop at a
#  P2 dir -d, wheels {2,3},           T2 turns then stop at b
#  P3 dir d,  wheel {3},              T3 turns then stop at c
#  P4 dir -d, no wheel (cam only), until a full turn (opening stroke)
# Each click enciphers one letter. Keystream variants: kd='wheels' (sum of wheel positions) or 'dial' (dial position).
# Unknown constant c0 (wheel start offsets) is free. K4 = clicks starting at offset o.
import numpy as np,itertools,random,time,sys
from stack import CT0,ANCH
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';AN=np.array(sorted(ANCH))
TURNS=[(3,2,1),(4,3,2),(2,1,0),(1,1,1),(0,0,0),(1,2,3)]
def stream(stops,turns,d,kd,start=0):
    p=start;w=[0,0,0];out=[]
    def click(dirn,eng):
        nonlocal p
        p=(p+dirn)%26
        for e in eng: w[e]=(w[e]+dirn)%26
        out.append((sum(w) if kd=='wheels' else p)%26)
    for ph,(stop,T,eng) in enumerate(zip(stops,turns,([0,1,2],[1,2],[2]))):
        dirn=d if ph%2==0 else -d
        n=26*T+((stop-p)*dirn)%26
        if n==0: n=26
        for _ in range(n): click(dirn,eng)
    for _ in range(26): click(-d if len(stops)%2==1 else d,[])
    s=np.array(out)
    while len(s)<97+len(out): s=np.concatenate([s,s])
    return s,len(out)
def Ks(ct):
    o={}
    for an,A in (('ST',ST),('KA',KA)):
        C=np.array([A.index(ct[i]) for i in AN]);P=np.array([A.index(ANCH[i]) for i in AN])
        o[(an,'vig')]=(C-P)%26;o[(an,'beau')]=(C+P)%26
    return o
def test_direct(ct,combos,ring_alpha_dep=True,starts=(0,)):
    K=Ks(ct);hits=[];n=0
    for (an,rule),k in K.items():
        A=ST if an=='ST' else KA
        for lab,stops_letters in combos:
            stops=[A.index(x) if isinstance(x,str) else x%26 for x in stops_letters]
            for turns in TURNS:
                for d in (1,-1):
                    for kd in ('wheels','dial'):
                        for st in starts:
                            s,L=stream(stops,turns,d,kd,st)
                            W=s[np.arange(L)[:,None]+AN[None,:]]
                            for sg in (1,-1):   # key = +stream or -stream
                                D=(k-sg*W)%26;ok=np.all(D==D[:,:1],axis=1);n+=L
                                for o in np.nonzero(ok)[0]: hits.append((lab,an,rule,turns,d,kd,st,sg,int(o)))
    return hits,n
if __name__=='__main__':
    t=time.time();SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    which=sys.argv[1]
    if which=='errors':
        C=[]
        for nm,trip in (('wrong A,U,Q (K3>K2>K1)',('A','U','Q')),('correct E,O,L (K3>K2>K1)',('E','O','L')),
                        ('wrong Q,U,A (K1>K2>K3)',('Q','U','A')),('correct L,O,E (K1>K2>K3)',('L','O','E')),
                        ('error positions K3,K2,K1 (10,114,56)',(10,114,56)),('error positions K1,K2,K3',(56,114,10)),
                        ('shift sizes K3>K2>K1',(None,))):
            if trip==(None,):
                for an,A in (('ST',ST),('KA',KA)):
                    C.append((nm+' '+an,tuple((A.index(b)-A.index(a))%26 for a,b in (('E','A'),('O','U'),('L','Q')))))
            else: C.append((nm,trip))
        combos=C;starts=range(26)
    elif which=='yar':
        combos=[('YAR',('Y','A','R')),('RAY',('R','A','Y')),('YAR+1',(25,1,18)),('RAY+1',(18,1,25))];starts=range(26)
    else:  # exhaustive, anchors decide
        combos=[(abc,abc) for abc in itertools.product(range(26),repeat=3)];starts=(0,)
    # control
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    s,L=stream([ST.index('A'),ST.index('U'),ST.index('Q')],(3,2,1),1,'wheels',4)
    ctl=''.join(ST[(ST.index(c)+s[i+40]+9)%26] for i,c in enumerate(pt))
    h,_=test_direct(ctl,[('ctrl',('A','U','Q'))],starts=(4,));print('control',[x for x in h if x[1]=='ST' and x[2]=='vig' and x[-1]==40][:1],flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        h,n=test_direct(ct,combos,starts=starts);print(which,nm,'window tests',n,'hits',len(h),h[:5],round(time.time()-t),flush=True)
