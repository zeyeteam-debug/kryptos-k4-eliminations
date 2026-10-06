from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
Wp=[i for i,c in enumerate(CT0) if c=='W'];starts=[0]+[w+1 for w in Wp]
keys=['QUA','QUAEU','LOE','LOEIN','LQOUEA','QLUOAE','YAR','RAY','YARQUA','QUAYAR','LOEYAR','DYAR','YARD','UQA','AUQ','EQUAL','QUAD']
# also shift-values from misspellings (L->Q etc) in ST and KA
best=[]
for kw in keys:
    for an,A in (('ST',ST),('KA',KA)):
        for rule in (0,1,2):
            for mode in ('global','restartW'):
                for off in range(len(kw)):
                    m=0
                    for i,p in ANCH.items():
                        j=i if mode=='global' else i-max(s for s in starts if s<=i)
                        k=A.index(kw[(j+off)%len(kw)]);c=A.index(CT0[i]);pp=A.index(p)
                        exp=[(pp+k)%26,(k-pp)%26,(pp-k)%26][rule]
                        m+= exp==c
                    best.append((m,kw,an,rule,mode,off))
best.sort(reverse=True);print('best matches out of 24:',best[:5])
# shift-sequence keys (numeric) from misspellings
for an,A in (('ST',ST),('KA',KA)):
    sh=[(A.index(b)-A.index(a))%26 for a,b in ('LQ','OU','EA','IE','NU')]
    print(an,'misspelling shifts',sh)
