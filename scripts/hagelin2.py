# Simplified Hagelin M-209 / C-38 tests on K4 cribs.
# key k_i = c0 + sum_j w_j * pin_j[(i+off_j) mod p_j]   (lug weights, no lug overlap), ST alphabet,
# Beaufort (k = c + p, the real M-209 rule up to the constant c0) and Vigenere (k = c - p) for safety.
import numpy as np,itertools,random,sys
from stack import CT0,ANCH,solve_mod
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
WHEELS=[26,25,23,21,19,17]
WL={26:ST,25:ST.replace('W',''),23:'ABCDEFGHIJKLMNOPQRSTUVX',21:ST[:21],19:ST[:19],17:ST[:17]}
idx=sorted(ANCH)
def kvals(ct,rule):
    return np.array([((ST.index(ct[i])+ST.index(ANCH[i])) if rule=='beau' else (ST.index(ct[i])-ST.index(ANCH[i])))%26 for i in idx])
# ---------- Part 1: free pins on two wheels, free weights & c0  (exact 2-SAT) ----------
def twosat(n,clauses):
    # literals: 2*v (true) / 2*v+1 (false)
    g=[[] for _ in range(2*n)]
    for a,b in clauses: g[a^1].append(b);g[b^1].append(a)
    idxc=[0];st=[];on=[False]*(2*n);low=[0]*(2*n);num=[-1]*(2*n);comp=[-1]*(2*n);c=[0]
    sys.setrecursionlimit(10000)
    def dfs(v):
        num[v]=low[v]=idxc[0];idxc[0]+=1;st.append(v);on[v]=True
        for w in g[v]:
            if num[w]<0: dfs(w);low[v]=min(low[v],low[w])
            elif on[w]: low[v]=min(low[v],num[w])
        if low[v]==num[v]:
            while True:
                w=st.pop();on[w]=False;comp[w]=c[0]
                if w==v: break
            c[0]+=1
    for v in range(2*n):
        if num[v]<0: dfs(v)
    return all(comp[2*v]!=comp[2*v+1] for v in range(n))
def pair_test(K,p1,p2):
    hits=[]
    W1,W2,C0=np.meshgrid(np.arange(1,26),np.arange(1,26),np.arange(26),indexing='ij')
    W1=W1.ravel();W2=W2.ravel();C0=C0.ravel()
    ok=np.ones(len(W1),bool)
    for v in K:
        d=(v-C0)%26;ok&=(d==0)|(d==W1)|(d==W2)|(d==(W1+W2)%26)
    for w1,w2,c0 in zip(W1[ok],W2[ok],C0[ok]):
        cl=[];bad=False
        for i,v in zip(idx,K):
            a=i%p1;b=p1+i%p2
            allowed=[(x,y) for x in (0,1) for y in (0,1) if (c0+w1*x+w2*y-v)%26==0]
            if not allowed: bad=True;break
            for x in (0,1):
                for y in (0,1):
                    if (x,y) not in allowed: cl.append((2*a+(1 if x==1 else 0),2*b+(1 if y==1 else 0)))
        if not bad and twosat(p1+p2,cl): hits.append((int(w1),int(w2),int(c0)))
    return hits
def part1(ct):
    out={}
    for p1,p2 in itertools.combinations(WHEELS,2):
        for rule in ('beau','vig'):
            out[(p1,p2,rule)]=len(pair_test(kvals(ct,rule),p1,p2))
    return out
# ---------- Part 2: pins fixed by a letter set (YAR, misspellings), weights linear, offsets enumerated ----------
PINSETS={'YAR':'YAR','QUAE':'QUAE','LOEI':'LOEI','QUAELOEI':'QUAELOEI','LQ+OU+EA+IE':'LQOUEAIE','KRYPTOS':'KRYPTOS'}
def part2(ct,pinsets=PINSETS,maxw=3,neg=True):
    hits=[];tested=0
    for nm,letters in pinsets.items():
        for comp in ((False,True) if neg else (False,)):
            pins={p:np.array([(c in letters)!=comp for c in WL[p]],int) for p in WHEELS}
            use=[p for p in WHEELS if 0<pins[p].sum()<p]
            for r in range(1,maxw+1):
                for sub in itertools.combinations(use,r):
                    if r==3 and not (19 in sub and 17 in sub): continue
                    for offs in itertools.product(*[range(p) for p in sub]):
                        for rule in ('beau','vig'):
                            K=kvals(ct,rule)
                            rows=[[int(pins[p][(i+o)%p]) for p,o in zip(sub,offs)]+[1] for i in idx]
                            o2,_=solve_mod(rows,list(K%2),2);o13,rk=solve_mod(rows,list(K%13),13);tested+=1
                            if o2 and o13: hits.append((24-rk,nm,comp,sub,offs,rule))
    return hits,tested
if __name__=='__main__':
    SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    # controls
    rng=random.Random(3);pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    b1=[rng.randrange(2) for _ in range(19)];b2=[rng.randrange(2) for _ in range(17)]
    ctl=''.join(ST[(5+9*b1[i%19]+14*b2[i%17]-ST.index(c))%26] for i,c in enumerate(pt))
    print('control part1 (19,17) beau w=9,14,c0=5 ->',(9,14,5) in pair_test(kvals(ctl,'beau'),19,17),flush=True)
    pins={p:[int(c in 'YAR') for c in WL[p]] for p in WHEELS}
    ctl2=''.join(ST[(3+7*pins[23][(i+4)%23]+11*pins[19][(i+10)%19]-ST.index(c))%26] for i,c in enumerate(pt))
    h,_=part2(ctl2,{'YAR':'YAR'},maxw=2,neg=False);print('control part2 ->',[x for x in h if x[3]==(23,19) and x[4]==(4,10)][:1],flush=True)
    for nm,ct in ():
        r=part1(ct);print('PART1',nm,{k:v for k,v in r.items() if v},'total sets',sum(r.values()),flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        h,t=part2(ct);h=[x for x in h if x[0]>=1]
        print('PART2',nm,'tested',t,'consistent',len(h),sorted(h,key=lambda x:-x[0])[:5],flush=True)
