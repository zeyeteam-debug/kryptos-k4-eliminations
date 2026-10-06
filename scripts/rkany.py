import random,itertools,numpy as np
from collections import Counter
exec(open('board.py').read().split('# K4 cell')[0])
from stack import CT0,ANCH
KA=KRY;ST=STD
AN=sorted(ANCH);PR=list(itertools.combinations(range(len(AN)),2))
CTpanel=''.join(ch for L in LINES for ch in L if ch.isalpha())
MORSE='VIRTUALLYINVISIBLEDIGETALINTERPRETATIUSHADOWFORCESLUCIDMEMORYTISYOURPOSITIONSOSRQ'
tab=''.join(KA[i:]+KA[:i] for i in range(26))
cities=('AZOREN KAPVERDE REYKJAVIK DUBLIN LONDON LISSABON MADEIRA BISSAU CASABLANCA CONAKRY DAKAR BAMAKO ACCRA '
 'AMSTERDAM BERLIN BRUSSEL BUDAPEST MADRID PARIS PRAG STOCKHOLM WARSCHAU OSLO KOPENHAGEN WIEN BERN PRESSBURG BELGRAD ROM TUNIS KINSHASA '
 'MOSKAU JEKATERINBURG BISCHKEK DUSCHANBE NEWDELHI KARACHI COLOMBO OMSK ALMATY TASCHKENT NOWOSIBIRSK RANGUN DHAKA KRASNOJARSK HANOI BANGKOK PHNOMPENH JAKARTA '
 'IRKUTSK ULANBATOR PEKING SHANGHAI MANILA PERTH HONGKONG KUALALUMPUR SINGAPUR JAKUTSK PJONGJANG TOKYO SEOUL CHABAROWSK WLADIWOSTOK SYDNEY CANBERRA MELBOURNE MAGADAN SACHALIN KAMTSCHATKA DATUMSGRENZE WELLINGTON').replace(' ','')
SRC={'K1P':K1P,'K2P':K2P,'K3P':K3P,'K123P':K1P+K2P+K3P,'CTpanel':CTpanel[:-97],'K4itself':CT0,'Morse':MORSE,'Tableau':tab,'Cities':cities,'KRYPTOSx':'KRYPTOS'*20,'K1-3+Morse':K1P+K2P+K3P+MORSE}
def keyvals(ct,A,rule):
    return [((A.index(ct[i])-A.index(ANCH[i])) if rule==0 else (A.index(ct[i])+A.index(ANCH[i])))%26 for i in AN]
def scan(ct):
    out=[]
    KV={(an,r):keyvals(ct,A,r) for an,A in (('ST',ST),('KA',KA)) for r in (0,1)}
    P=[ANCH[i] for i in AN];C=[ct[i] for i in AN]
    for nm,s in SRC.items():
        for dn,t in (('fwd',s),('rev',s[::-1])):
            L=len(t);ext=t*(97//L+3)
            for off in range(L):
                kl=[ext[off+i] for i in AN]
                samek=[(a,b) for a,b in PR if kl[a]==kl[b]]
                if not samek: continue
                # (b) alphabet ST/KA, each key letter -> arbitrary shift: same key letter => same shift
                for key,K in KV.items():
                    if all(K[a]==K[b] for a,b in samek): out.append(('shift',len(samek),nm,dn,off)+key)
                # (a) fully arbitrary tableau: same key letter: same p <=> same c
                strong=sum(P[a]==P[b] for a,b in samek)
                if all((P[a]==P[b])==(C[a]==C[b]) for a,b in samek) and strong>=1: out.append(('anytab',strong,len(samek),nm,dn,off))
    return out
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list((K3P+K1P)[10:107]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    key=(K2P*2)[77:174];perm=random.Random(3).sample(ST,26)   # random mixed alphabet tableau
    ctrl=''.join(perm[(perm.index(p)+ST.index(k))%26] for p,k in zip(pt,key))
    h=scan(ctrl);print('control (K2P off 77, random alphabet): found',[x for x in h if x[2:5]==('K2P','fwd',77) or x[3:6]==('K2P','fwd',77)][:3])
    print('control hit counts',Counter(x[0] for x in h))
    nulls=Counter()
    for sd in range(20):
        l=list(CT0);random.Random(sd).shuffle(l)
        for x in scan(''.join(l)):
            nulls[(x[0],x[1])]+=1
    print('null (20 shuffles) counts by (type,strength):',sorted(nulls.items()))
    hk=scan(CT0);print('K4 counts by (type,strength):',sorted(Counter((x[0],x[1]) for x in hk).items()))
    print('K4 strongest:',sorted(hk,key=lambda x:-x[1])[:8])
