import numpy as np
import dictsearch as D
from stack import ANCH
from sculptalpha import ALS
AL=ALS['PALIMPSEST/ST']
pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
for i,c in ANCH.items(): pt[i]=c
x=[AL[(AL.index(c)+AL.index('SHADOW'[i%6]))%26] for i,c in enumerate(pt)]
m=D.perm_from_read(D.spiral_order(97,9,'BL',1),97);y=['?']*97
for i in range(97): y[m[i]]=x[i]
TB=D.table(D.WORDS,ANCH,AL);h=D.run(''.join(y),ANCH,[('sp',m)],[],AL,TB)
print('control (PALIMPSEST alphabet, Q(SHADOW), spiral 9 BL):',[r for r in h if r[0]=='SHADOW'][:1])
