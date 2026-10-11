import numpy as np,random,string
from stack import CT0,ANCH
from gen import run
from cribwords import BW,free_slots
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
def gt(ct,words,P=19,S=2):
    res={}
    for w in words:
        n=h=0;strong=[]
        for s in free_slots(len(w)):
            an=dict(ANCH);an.update({s+j:c for j,c in enumerate(w)})
            H=[x for x in run(ct,an,np.arange(97,dtype=np.int16)[None,:],[(w,s)],Pmax=P,anyalpha=True,ANYSTRONG=S,RTHR=99) if x[0]=='any' and x[2]==P]
            n+=1
            if H: h+=1;strong.append((s,max(x[3] for x in H)))
        res[w]=(n,h,strong)
    return res
rng=random.Random(11)
pt=list('LENINPLATZTHEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNK'[:97])
for i,c in ANCH.items(): pt[i]=c
alph=[''.join(rng.sample(ST,26)) for _ in range(19)]
ctl=''.join(alph[i%19][ST.index(c)] for i,c in enumerate(pt))
print('control',gt(ctl,['LENINPLATZ'])['LENINPLATZ'])
r=gt(CT0,BW);N=sum(v[0] for v in r.values());H=sum(v[1] for v in r.values())
rr=random.Random(5);RW=[''.join(rr.choice(string.ascii_uppercase) for _ in range(len(w))) for w in BW for _ in range(3)]
r0=gt(CT0,RW);N0=sum(v[0] for v in r0.values());H0=sum(v[1] for v in r0.values())
print('K4 Berlin words: placements',N,'consistent',H,'rate',round(H/N,3),'| random words',N0,H0,'rate',round(H0/N0,3))
# rank Berlin words by strength (number of coincidences) beyond random
best=sorted([(st,w,s) for w,v in r.items() for s,st in v[2]],reverse=True)[:10]
print('strongest Berlin placements',best)
print('strongest random',sorted([(st,w,s) for w,v in r0.items() for s,st in v[2]],reverse=True)[:5])
