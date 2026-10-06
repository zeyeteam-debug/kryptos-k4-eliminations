import itertools,sys
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
CT0="OBKRUOXOGHULBSOLIFBBWFLRVQQPRNGKSSOTWTQSJQSSEKZZWATJKLUDIAWINFBNYPVTTMZFPKWGDKZXTJCDIGKUHUAUEKCAR"
ANCH={}
for s,w in [(21,'EASTNORTHEAST'),(63,'BERLINCLOCK')]:
    for j,c in enumerate(w): ANCH[s+j]=c
def solve_mod(rows,rhs,m):
    # gaussian elimination over prime field m; returns (consistent, rank)
    A=[r[:]+[b] for r,b in zip(rows,rhs)];n=len(rows[0]);rank=0;r=0
    for col in range(n):
        piv=None
        for i in range(r,len(A)):
            if A[i][col]%m: piv=i;break
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r];inv=pow(A[r][col],m-2,m)
        A[r]=[(x*inv)%m for x in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][col]%m:
                f=A[i][col];A[i]=[(x-f*y)%m for x,y in zip(A[i],A[r])]
        r+=1
    ok=all(any(x%m for x in row[:-1]) or row[-1]%m==0 for row in A)
    return ok,r
def test(CT,periods,AL,model):
    idx=sorted(ANCH);offs=[];o=0
    for p in periods: offs.append(o);o+=p
    rows=[];rhs=[]
    for i in idx:
        row=[0]*o
        for p,off in zip(periods,offs): row[off+i%p]+=1
        c=AL.index(CT[i]);pp=AL.index(ANCH[i])
        d=(c-pp)%26 if model=='vig' else (c+pp)%26
        rows.append(row);rhs.append(d)
    ok2,r2=solve_mod(rows,rhs,2);ok13,r13=solve_mod(rows,rhs,13)
    red=len(idx)-r13   # redundant equations (mod 13 rank)
    return ok2 and ok13,red
def campaign(CT):
    hits=[];expected=0.0;tested=0
    sets=[(p,) for p in range(1,25)]+[t for t in itertools.combinations(range(2,16),2)]+[t for t in itertools.combinations(range(2,11),3)]+[t for t in itertools.combinations(range(2,9),4)]
    for per in sets:
        for AL,an in ((KA,'KA'),(ST,'ST')):
            for model in ('vig','beau'):
                ok,red=test(CT,per,AL,model);tested+=1
                if red>=1: expected+=26.0**(-red)
                if ok and red>=3: hits.append((red,per,an,model))
    return hits,expected,tested
if __name__=='__main__':
    h,e,t=campaign(CT0)
    print('tested',t,'expected false positives (red>=1)',round(e,3))
    print('consistent with >=3 redundant checks:',sorted(h,reverse=True)[:20])
    # positive control: plant sum of keys OCHRE + VERDIGRIS
    import random
    pt=list('X'*97);src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK'
    ct=''.join(KA[(KA.index(c)+KA.index('OCHRE'[i%5])+KA.index('VERDIGRIS'[i%9]))%26] for i,c in enumerate(pt))
    print('control (5,9) KA vig ->',test(ct,(5,9),KA,'vig'))
    ct3=''.join(KA[(KA.index(c)+KA.index('METE'[i%4])+KA.index('METER'[i%5])+KA.index('METIER'[i%6]))%26] for i,c in enumerate(pt))
    print('control (4,5,6) KA vig ->',test(ct3,(4,5,6),KA,'vig'))
