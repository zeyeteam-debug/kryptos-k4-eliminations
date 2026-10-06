import random
from decimal import Decimal,getcontext
from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
AN=sorted(ANCH);D=600
getcontext().prec=D+10
def pi_digits(n):
    def arctan_inv(x,u):
        s=t=u//x;x2=x*x;k=1;sign=-1
        while t: t//=x2;s+=sign*(t//(2*k+1));sign=-sign;k+=1
        return s
    u=10**(n+10);p=4*(4*arctan_inv(5,u)-arctan_inv(239,u));return str(p)[:n]
def dec_digits(x,n): return str(x).replace('.','')[:n]
e=sum(Decimal(1)/Decimal(__import__('math').factorial(k)) for k in range(400))
consts={'pi':pi_digits(D),'e':dec_digits(e,D),'sqrt2':dec_digits(Decimal(2).sqrt(),D),'sqrt3':dec_digits(Decimal(3).sqrt(),D),
 'phi':dec_digits((1+Decimal(5).sqrt())/2,D)}
# digit streams -> numeric key sequences
seqs={}
for nm,ds in consts.items():
    seqs[nm+'_1digit']=[int(c) for c in ds]
    seqs[nm+'_2digit']=[int(ds[i:i+2]) for i in range(0,len(ds)-1,2)]
    seqs[nm+'_2digit_sliding']=[int(ds[i:i+2]) for i in range(len(ds)-1)]
dates={'wall_11091989':'11091989','wall_09111989':'09111989','dedic_03111990':'03111990','dedic_11031990':'11031990','coords':'385765778844','coords2':'3857657784','ene':'675','clock_19891109':'19891109'}
for nm,ds in dates.items(): seqs['date_'+nm]=[int(c) for c in ds]*40
n=300
seqs['squares']=[i*i for i in range(n)];seqs['cubes']=[i**3 for i in range(n)];seqs['triangular']=[i*(i+1)//2 for i in range(n)]
P=[];k=2
while len(P)<n:
    if all(k%p for p in P if p*p<=k): P.append(k)
    k+=1
seqs['primes']=P;seqs['prime_gaps']=[P[i+1]-P[i] for i in range(n-1)]
for b in (2,3,5,7):seqs[f'pow{b}']=[pow(b,i,26*1000) for i in range(n)]
L=[2,1]
for i in range(n): L.append(L[-1]+L[-2])
seqs['lucas']=L
def scan(ct):
    best=[];hits=[]
    for nm,s in seqs.items():
        for mul in (1,-1):
            for an,A in (('ST',ST),('KA',KA)):
                for rule in (0,1):
                    for off in range(0,max(1,len(s)-97)):
                        m=0
                        for i in AN:
                            k=(mul*s[off+i])%26;p=A.index(ANCH[i])
                            exp=(p+k)%26 if rule==0 else (k-p)%26
                            m+=A[exp]==ct[i]
                        best.append(m)
                        if m==24: hits.append((nm,mul,an,rule,off))
    return hits,max(best),len(best)
if __name__=='__main__':
    # control: pi 2-digit key from offset 5, KA add
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    s=seqs['pi_2digit']
    ctrl=''.join(KA[(KA.index(p)+s[5+i])%26] for i,p in enumerate(pt))
    print('control:',scan(ctrl)[0][:2])
    h,b,n=scan(CT0);print('K4 full matches:',h,'| best',b,'/24 | tested',n)
    l=list(CT0);random.Random(2).shuffle(l);print('null best',scan(''.join(l))[1])
