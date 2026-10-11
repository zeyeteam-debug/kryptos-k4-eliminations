import numpy as np
from stack import ANCH
from allalpha import *
A=ALS['PALIMPSEST/ST']
pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
for i,c in ANCH.items(): pt[i]=c
# A: clock key in PALIMPSEST alphabet
ct=''.join(A[(A.index(c)+5*i+11*((i+3)//12)+7)%26] for i,c in enumerate(pt))
labs,K=Kmat(ct);h=linear_test(K,[m for m in linear_models() if m[0]==('clock',12,3)]);print('A ctrl',[(r,l,labs[j]) for r,l,j in h][:2])
# D: running key from Carter text in PALIMPSEST alphabet (Quagmire III), offset 5000
c=''.join(ch for ch in open('texts/carter_vol1.txt',errors='ignore').read().upper() if 'A'<=ch<='Z')
ct=''.join(A[(A.index(p)+A.index(c[5000+i]))%26] for i,p in enumerate(pt))
print('D ctrl',[x for x in rk_test(ct,{'carter':c}) if x[0]=='PALIMPSEST/ST'][:2])
# F: ct-autokey L=7 in PALIMPSEST alphabet
ctl=[];prim='SHADOWS'
for i,p in enumerate(pt):
    k=prim[i] if i<7 else ctl[i-7];ctl.append(A[(A.index(p)+A.index(k))%26])
print('F ctrl',[x for x in autokey_test(''.join(ctl)) if x[0]=='PALIMPSEST/ST'][:2])
