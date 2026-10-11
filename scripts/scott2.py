# How surprising is "BERLIN somewhere" under Q3(KRYPTOS, 11-letter window of the sculpture ciphertext) + affine position map?
import numpy as np
from stack import CT0
from board import LINES
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';N=97
allct=''.join(LINES).replace('?','')
W=sorted({allct[i:i+11] for i in range(len(allct)-10)})
tgt=np.array([KA.index(c) for c in 'BERLIN']);C=np.array([KA.index(c) for c in CT0])
anyhit=0;poshit=0;keys_any=0
for w in W:
    k=np.array([KA.index(c) for c in w])
    for sign in (1,-1):
        s=(C-sign*k[np.arange(N)%11])%26
        hit=False
        for b in range(1,N):
            idx=(np.arange(N)[:,None]+b*np.arange(6)[None,:])%N     # t + b*j
            m=np.all(s[idx]==tgt,axis=1);c=int(m.sum());anyhit+=c;hit|=c>0
        keys_any+=hit
print('keyword windows',len(W),'x 2 signs; affine maps 97x96 per direction')
print('total (keyword, b, start) giving BERLIN somewhere:',anyhit,'; keywords with at least one such map:',keys_any)
