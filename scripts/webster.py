import random
from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
def keyed(w):
    s=''
    for c in w+ST:
        if c not in s: s+=c
    return s
Wp=[i for i,c in enumerate(CT0) if c=='W'];starts=[0]+[w+1 for w in Wp]
AN=sorted(ANCH)
names=['WILLIAMWEBSTER','WILLIAMHWEBSTER','WILLIAMHEDGCOCKWEBSTER','WEBSTER','WILLIAM','WW','BILLWEBSTER','JUDGEWEBSTER','WEBSTERWILLIAM']
def score(ct):
    best=[]
    for nm in names:
        alphs={'ST':ST,'KA':KA,'keyed':keyed(nm),'keyedKRYPTOS':keyed(nm+'KRYPTOS')}
        for pa,PA in alphs.items():
            for ca,CA in alphs.items():
                for rule in range(3):
                    for mode in ('global','restartW'):
                        for off in range(len(nm)):
                            m=0
                            for i in AN:
                                j=i if mode=='global' else i-max(s for s in starts if s<=i)
                                k=CA.index(nm[(j+off)%len(nm)]);p=PA.index(ANCH[i]);c=CA.index(CT0[i]) if False else None
                                exp=[(p+k)%26,(k-p)%26,(p-k)%26][rule]
                                m+= CA[exp]==ct[i]
                            best.append((m,nm,pa,ca,rule,mode,off))
    best.sort(reverse=True);return best
b=score(CT0);print('K4 best:',b[:4])
l=list(CT0);random.Random(1).shuffle(l);print('null best:',score(''.join(l))[0])
print('combinations tested:',len(b))
