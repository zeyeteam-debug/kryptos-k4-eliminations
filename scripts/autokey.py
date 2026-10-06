import random
from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
AN=sorted(ANCH)
def keyat(ct,A,rule):
    k={}
    for i in AN:
        c=A.index(ct[i]);p=A.index(ANCH[i])
        k[i]=(c-p)%26 if rule==0 else (c+p)%26 if rule==1 else (p-c)%26
    return k
def tests(ct):
    hits=[]
    for an,A in (('ST',ST),('KA',KA)):
        for rule in (0,1,2):
            k=keyat(ct,A,rule)
            # (a) ciphertext autokey, lag L, shift s
            for L in range(1,97):
                d=[(k[i]-A.index(ct[i-L]))%26 for i in AN if i>=L]
                if len(d)>=20 and len(set(d))==1: hits.append(('CTautokey',an,rule,L,len(d)))
                d=[(k[i]+A.index(ct[i-L]))%26 for i in AN if i>=L]   # reversed-sign key
                if len(d)>=20 and len(set(d))==1: hits.append(('CTautokey-',an,rule,L,len(d)))
            # (b) plaintext autokey, lag L
            for L in range(1,97):
                d=[(k[i]-A.index(ANCH[i-L]))%26 for i in AN if i-L in ANCH]
                if len(d)>=5 and len(set(d))==1: hits.append(('PTautokey',an,rule,L,len(d)))
                d=[(k[i]+A.index(ANCH[i-L]))%26 for i in AN if i-L in ANCH]
                if len(d)>=5 and len(set(d))==1: hits.append(('PTautokey-',an,rule,L,len(d)))
            # (c) Fibonacci / lagged recurrence key  k_i = x*k_{i-a} + y*k_{i-b} + s
            for a in range(1,12):
                for b in range(a+1,13):
                    for x in (1,25):
                        for y in (1,25):
                            d=[(k[i]-x*k[i-a]-y*k[i-b])%26 for i in AN if i-a in k and i-b in k]
                            if len(d)>=4 and len(set(d))==1: hits.append(('recur',an,rule,a,b,x,y,len(d)))
            # (d) Gromark: key digits 0-9 with mixed CT alphabet A, plain ST
            g=[(A.index(ct[i])-ST.index(ANCH[i]))%26 for i in AN]
            if max(g)<=9: hits.append(('gromark',an))
    return hits
def enc_ctauto(pt,A,L,primer,rule):
    ct=[]
    for i,ch in enumerate(pt):
        k=A.index(primer[i]) if i<L else A.index(ct[i-L])
        p=A.index(ch);c=(p+k)%26 if rule==0 else (k-p)%26
        ct.append(A[c])
    return ''.join(ct)
def enc_ptauto(pt,A,L,primer,rule):
    ct=[]
    for i,ch in enumerate(pt):
        k=A.index(primer[i]) if i<L else A.index(pt[i-L])
        p=A.index(ch);c=(p+k)%26 if rule==0 else (k-p)%26
        ct.append(A[c])
    return ''.join(ct)
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    print('control CT-autokey KA L=7:',[h for h in tests(enc_ctauto(pt,KA,7,'KRYPTOS',0)) if h[0].startswith('CT')][:2])
    print('control PT-autokey ST L=4 beaufort:',[h for h in tests(enc_ptauto(pt,ST,4,'CLOK',1)) if h[0].startswith('PT')][:2])
    print('K4:',tests(CT0))
    random.seed(1);fp=0
    for _ in range(300):
        s=list(CT0);random.shuffle(s);fp+=bool(tests(''.join(s)))
    print('null shuffles with any hit: %d/300'%fp)
