# Left panel (K1-K2 ciphertext / plaintext) laid as a mask with its top-left at Y, A or R: coincidence holes
exec(open('handmask.py').read().split("print('selections'")[0].replace("print(","(lambda *a,**k:None)("))
import random
SEL2={}
for an in ('Y','A','R','panel2_corner'):
    ar,ac=ANCHORS[an]
    for bn,B in (('ct',LINES),('pt',PB)):
        for mir in (False,True):
            for rel in ('eq','neq_none'):
                S=set()
                for i,(Lr,c) in enumerate(cells):
                    r,cc=Lr-ar,c-ac
                    if 0<=r<14:
                        row=B[r][::-1] if mir else B[r]
                        if 0<=cc<len(row) and ((rel=='eq' and row[cc]==CT0[i]) or (rel=='neq_none' and row[cc]=='?')): S.add(i)
                if 0<len(S)<N: SEL2[(an,bn,mir,rel)]=S
print({k:sorted(v) for k,v in SEL2.items()})
P2=[];L2=[]
for k,S in SEL2.items():
    a=[i for i in range(N) if i in S];b=[i for i in range(N) if i not in S]
    for nm,o in (('first',a+b),('last',b+a)):
        o=np.array(o);P2+=[o,np.argsort(o)];L2+=[(k,nm,'a'),(k,nm,'b')]
P2=np.array(P2,dtype=np.int16)
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
for nm,x in (('K4',CT0),('null',SH)):
    print(nm,len(run(x,ANCH,P2,L2,Pmax=24)),qii(x,P2,24))
