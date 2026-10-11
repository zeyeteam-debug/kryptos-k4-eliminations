# Arbitrary nulls anywhere (any number, any positions outside the cribs) + periodic key on the reduced text.
# Only two numbers matter for the cribs: a = nulls before pos 21, b = nulls between 34 and 62.
import random
from stack import CT0,ANCH,KA,ST
def ok(ct,a,b,p,mode):
    cl={}
    for i,P in ANCH.items():
        r=i-a if i<34 else i-a-b
        k=r%p;C=ct[i]
        if mode=='any':
            d=cl.setdefault(k,{});e=cl.setdefault(('inv',k),{})
            if d.get(P,C)!=C or e.get(C,P)!=P: return False
            d[P]=C;e[C]=P
        else:
            A=KA if mode[0]=='K' else ST;s=1 if mode[1]=='v' else -1
            v=(A.index(C)-s*A.index(P))%26
            if cl.setdefault(k,v)!=v: return False
    return True
def red(a,b,p): # need >=6 repeated-class constraints to be informative
    rs=[(i-a if i<34 else i-a-b)%p for i in ANCH];return len(rs)-len(set(rs))
def scan(ct):
    H=[]
    for a in range(0,22):
        for b in range(0,30):
            for p in range(1,27):
                if red(a,b,p)<6: continue
                for m in ('Kv','Kb','Sv','Sb','any'):
                    if ok(ct,a,b,p,m): H.append((a,b,p,m))
    return H
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
# control: plaintext with 5 nulls before 21 and 7 between cribs, KA vig period 9
pt=[c for c in 'THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTO']
full=['?']*97
for i,c in ANCH.items(): full[i]=c
real=[i for i in range(97) if not( i in (2,5,9,12,17) or i in (36,40,44,48,52,56,60))]
ct=['?']*97;j=0
for r,i in enumerate(real):
    P=full[i] if full[i]!='?' else pt[r%len(pt)]
    ct[i]=KA[(KA.index(P)+KA.index('ABSCISSAX'[r%9]))%26]
for i in range(97):
    if ct[i]=='?': ct[i]=random.choice(ST)
print('control',[h for h in scan(''.join(ct)) if h[:3]==(5,7,9)])
for nm,x in (('K4',CT0),('null',SH)):
    H=scan(x);print(nm,len(H),H[:10])
