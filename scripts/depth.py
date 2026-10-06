from stack import CT0,ANCH
Wp=[i for i,c in enumerate(CT0) if c=='W'];starts=[0]+[w+1 for w in Wp];ends=Wp+[97]
segs=[(s,e) for s,e in zip(starts,ends)]
# alphabet per local offset j: ct->pt known pairs
A={}
for i,p in ANCH.items():
    s=max(x for x in starts if x<=i);j=i-s;A.setdefault(j,{})[CT0[i]]=p
for s,e in segs:
    ct=CT0[s:e];pt=''.join(A.get(j,{}).get(c,'.') for j,c in enumerate(ct))
    print(f'{s:2d}-{e-1:2d}  CT {ct}\n        PT {pt}')
