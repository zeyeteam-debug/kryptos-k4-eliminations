# Punched card (Hollerith, 12 rows x 80 cols) as a grille. Word punched on the card -> hole mask; K4 laid on the card grid;
# read letters under holes first, then the rest (and inverse). Then any periodic key <=24 (ST/KA/any alphabet) + clock/gear models.
import numpy as np,random,time
from stack import CT0,ANCH
from gen import run
ROWS=['12','11','0','1','2','3','4','5','6','7','8','9']
def holl(ch):
    if ch.isdigit(): return {ch}
    i=ord(ch)-65
    if i<9: return {'12',str(i+1)}
    if i<18: return {'11',str(i-8)}
    return {'0',str(i-16)}
def mask(word):
    return np.array([[1 if r in holl(c) else 0 for c in word] for r in ROWS])   # 12 x len
def perms_for(word):
    m=mask(word);L=m.shape[1];out=[]
    for layout in ('colmajor','rowmajor'):
        if layout=='colmajor':   # K4 down the card columns (12 letters per column), mask tiled across columns
            cells=[(r,c) for c in range(9) for r in range(12)][:97]
        else:                    # K4 across rows of width 9 (12 rows x 9 cols >= 97)
            cells=[(r,c) for r in range(11) for c in range(9)][:97]
        for shift in range(L):
            h=[i for i,(r,c) in enumerate(cells) if m[r,(c+shift)%L]];o=[i for i in range(97) if i not in set(h)]
            if 0<len(h)<97:
                order=np.array(h+o);out.append(((word,layout,shift,'a'),order));out.append(((word,layout,shift,'b'),np.argsort(order)))
    return out
if __name__=='__main__':
    from sculptalpha import words as SW
    from clockwords import W as CW
    from greek import G
    words=sorted(set(SW)|set(CW)|set(G)|{'KRYPTOS','PALIMPSEST','ABSCISSA','BERLINCLOCK','EASTNORTHEAST','CIA','LANGLEY','SCHEIDT','SANBORN','NOVEMBER1990','1990','3NOV1990','11091989'})
    words=[w for w in words if w.isalnum()]
    P=[];L=[]
    for w in words:
        for lab,o in perms_for(w): P.append(o);L.append(lab)
    P=np.array(P,dtype=np.int16);print('words',len(words),'grille perms',len(P),flush=True)
    import mech
    from mech import check,special
    keep_r={12,24,60}|{x for t in special.values() for x in t if x>=2}
    sel=[(l,A) for (l,A) in mech.ALL if (l[0]=='clock' and l[1] in keep_r) or l[0] in('clock3','gear')]
    mech.build(sel);s2=[m for m,r in zip(sel,mech.RED) if r>=6];mech.build(s2)
    # control: KRYPTOS card colmajor shift 0, then KA vig period 7
    KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';lab,o=perms_for('KRYPTOS')[0]
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    x=[KA[(KA.index(c)+KA.index('SHADOWS'[i%7]))%26] for i,c in enumerate(pt)];m=np.argsort(o);y=['?']*97
    for i in range(97): y[m[i]]=x[i]
    print('control',[h for h in run(''.join(y),ANCH,np.array([m],dtype=np.int16),[lab],Pmax=12) if h[0]=='key'][:1],flush=True)
    SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    for nm,ct in (('K4',CT0),('null',SH)):
        H=run(ct,ANCH,P,L,Pmax=24);h2,e=check(ct,P,L,minred=6)
        print(nm,'periodic/any-alpha hits',len(H),H[:3],'| clock/gear',len(h2),flush=True)
