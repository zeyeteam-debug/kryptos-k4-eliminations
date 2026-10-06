import random
from stack import solve_mod,ANCH,CT0
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
idx=sorted(ANCH)
def check(ct,rowf,nun,AL,rule):
    rows=[];rhs=[]
    for i in idx:
        rows.append(rowf(i));c=AL.index(ct[i]);p=AL.index(ANCH[i])
        rhs.append((c-p)%26 if rule==0 else (c+p)%26)
    ok2,r2=solve_mod(rows,rhs,2);ok13,r13=solve_mod(rows,rhs,13)
    return ok2 and ok13, len(idx)-max(r2,r13)
def models():
    # M1 odometer: minute wheel period p, hour wheel q positions advancing every p letters
    for p in range(2,25):
        for q in range(2,25):
            if p+q>20: continue
            def f(i,p=p,q=q):
                r=[0]*(p+q);r[i%p]+=1;r[p+(i//p)%q]+=1;return r
            yield ('odometer',p,q),f,p+q
    # M2 Alberti-style stepping: periodic key p plus shift d every N letters
    for p in range(1,19):
        for N in range(2,49):
            def f(i,p=p,N=N):
                r=[0]*(p+1);r[i%p]+=1;r[p]=i//N;return r
            yield ('alberti',p,N),f,p+1
    # M3 24-hour clock: key = a[i mod p] + hour(i) where hour advances every N letters mod 24 (values not mod 26!)
def run(ct,thr):
    hits=[]
    for name,f,n in models():
        for an,AL in (('ST',ST),('KA',KA)):
            for rule in (0,1):
                ok,red=check(ct,f,n,AL,rule)
                if ok and red>=thr: hits.append((red,name,an,rule))
    return sorted(hits,reverse=True)
def plant(pt,AL,key):
    return ''.join(AL[(AL.index(ch)+key(i))%26] for i,ch in enumerate(pt))
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    a=[KA.index(x) for x in 'CLOCK'];b=[KA.index(x) for x in 'BERLIN']
    c1=plant(pt,KA,lambda i:a[i%5]+b[(i//5)%6])
    print('ctrl odometer(5,6):',[h for h in run(c1,3)][:3])
    c2=plant(pt,ST,lambda i:[3,1,4][i%3]+7*(i//12))
    print('ctrl alberti(3,N=12):',[h for h in run(c2,3)][:3])
    random.seed(2);nulls=[]
    for _ in range(30):
        s=list(CT0);random.shuffle(s);h=run(''.join(s),3);nulls.append(h[0][0] if h else 0)
    print('null best redundancy (30 shuffles):',sorted(nulls))
    print('K4:',run(CT0,3)[:10])
