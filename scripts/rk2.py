from rkany import *
P=[ANCH[i] for i in AN];C=[CT0[i] for i in AN]
print('K4 anchor pairs with same PT and same CT:',[(AN[a]+1,AN[b]+1,P[a],C[a]) for a,b in PR if P[a]==P[b] and C[a]==C[b]])
cnt=[]
for sd in range(200):
    l=list(CT0);random.Random(1000+sd).shuffle(l);cnt.append(sum(1 for x in scan(''.join(l)) if x[0]=='anytab'))
import numpy as np;cnt=np.array(cnt)
print('null anytab hits: mean',cnt.mean(),'95%',np.percentile(cnt,95),'max',cnt.max(),'share >=143:',(cnt>=143).mean())
# shift-type control: KA tableau but key letters mapped to arbitrary shifts
src=(K3P+K1P)[10:107];pt=list(src);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
key=(K1P+K2P+K3P)[200:297];sh={c:random.Random(7).randrange(26)*0+((ord(c)*7)%26) for c in ST}
ctrl=''.join(KA[(KA.index(p)+sh[k])%26] for p,k in zip(pt,key))
print('shift control:',[x for x in scan(ctrl) if x[0]=='shift'][:3])
