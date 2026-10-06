import random
from stack import CT0,ANCH
def check(ct,sep):
    Wp=[i for i,c in enumerate(ct) if c==sep];starts=[0]+[w+1 for w in Wp]
    if any(ct[i]==sep for i in ANCH): return None
    loc={i:i-max(s for s in starts if s<=i) for i in ANCH}
    pairs=[(a,b) for a in range(21,34) for b in range(63,74) if loc[a]==loc[b]]
    if not pairs: return None
    bad=sum((ANCH[a]==ANCH[b])!=(ct[a]==ct[b]) for a,b in pairs);strong=sum(ANCH[a]==ANCH[b] for a,b in pairs)
    return bad,strong,len(pairs)
print('K4 per separator letter:')
for L in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
    r=check(CT0,L)
    if r: print(' ',L,'contradictions=%d strong=%d overlaps=%d'%r)
# null: shuffle non-W letters keeping W positions
rng=random.Random(0);n=0;ok=0
nonW=[i for i,c in enumerate(CT0) if c!='W'];letters=[CT0[i] for i in nonW]
for _ in range(20000):
    rng.shuffle(letters);s=list(CT0)
    for i,c in zip(nonW,letters): s[i]=c
    r=check(''.join(s),'W')
    if r: n+=1;ok+= (r[0]==0 and r[1]>=1)
print('null (W kept, rest shuffled): pass with 0 contradictions & >=1 strong: %.3f'%(ok/n))
