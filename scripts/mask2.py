import mask,random
from stack import CT0
KA=mask.KA
# sliding mask over tableau: equivalent to tableau whose rows shift by alpha per row
G={}
for a in range(26):
    T=[''.join(KA[(a*r+j)%26] for j in range(32)) for r in range(28)]
    G[f'slide{a}']=T
mask.GRIDS=G
h=mask.scan(CT0);print('K4 sliding-mask hits',len(h),h[:5])
l=list(CT0);random.Random(9).shuffle(l);print('null',len(mask.scan(''.join(l))))
