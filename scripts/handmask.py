# Hand-made masks derivable "from home": raised YAR as the mask's top-left corner, Morse-E masks,
# left-panel features, tableau coincidences, text-derived masks, Josephus skip readings.
import numpy as np,random,time,itertools
exec(open('board.py').read().split('best=[];trials=0')[0])
from stack import CT0,ANCH,KA,ST
from gen import run
from qii_more import scan as qii
from k3rot import compose
import os
TP='tableau.txt' if os.path.exists('tableau.txt') else 'repos/keithadler_kryptos/data/tableau.txt'  # tableau as carved (Wikipedia, CC BY-SA)
TAB=[l.strip() for l in open(TP) if l.strip()]
N=97;CR=set(ANCH)
cellidx={c:i for i,c in enumerate(cells)}
# letter cells of the whole cipher side in reading order (skipping ?)
LC=[(L,c) for L,l in enumerate(LINES) for c,ch in enumerate(l) if ch!='?']
# misspelling cells
q=K1P.index('IQLUSION')+1; k2=K2P.index('UNDERGRUUND')+7
K3C=''.join(LINES)[''.join(LINES).index('ENDYAH'):''.join(LINES).index('?OBKR')].replace('?','')
idx=compose(336,42,True,14,True); p3=K3P.index('DESPARATLY')+4; n3=int(np.nonzero(idx==p3)[0][0])
assert ''.join(np.array(list(K3P))[idx])==K3C
off3=len(K1P)+len(K2P)
MIS={'IQLUSION_Q':LC[q],'UNDERGRUUND_U':LC[len(K1P)+k2],'DESPARATLY_A':LC[off3+n3]}
QM=[(L,c) for L,l in enumerate(LINES) for c,ch in enumerate(l) if ch=='?']
ANCHORS={'Y':(14,3),'A':(14,4),'R':(14,6),'panel2_corner':(14,0),'sheet_corner':(0,0)}
for i,qq in enumerate(QM): ANCHORS[f'?{i}']=qq
ANCHORS.update(MIS)
assert LINES[14][3]=='Y' and LINES[14][4]=='A' and LINES[14][6]=='R'
# ---- templates: list of relative (dr,dc) holes
def morse_rows(split):
    lines=[[('E',2),'VIRTUALLY',('E',1),'|',('E',6),'INVISIBLE'],['DIGETAL',('E',3),'|','INTERPRETATIT'],
           [('E',2),'SHADOW',('E',2),'|','FORCES',('E',5)],['LUCID',('E',3),'|','MEMORY',('E',1)],['TISYOUR','|','POSITION',('E',1)],['SOS'],['RQ']]
    rows=[]
    for ln in lines:
        segs=[[]]
        for t in ln:
            if t=='|':
                if split: segs.append([])
                continue
            segs[-1]+=([1]*t[1]) if isinstance(t,tuple) else [0]*len(t)
        rows+=segs
    return rows
TEMPL={}
for split in (False,True):
    R=morse_rows(split)
    for inv in (False,True):
        TEMPL[f'morseE_split{split}_inv{inv}']=[(r,c) for r,row in enumerate(R) for c,v in enumerate(row) if v!=inv]
TEMPL['left_?']=[(L,c) for L,c in QM if L<14]
TEMPL['left_misspell']=[MIS['IQLUSION_Q'],MIS['UNDERGRUUND_U']]
TEMPL['all_misspell_rel_panel']=[(L%14,c) for L,c in MIS.values()]
TEMPL['YAR_shape']=[(0,0),(0,1),(0,3)]
for h in range(1,9):
    for w in range(1,9):
        if h*w>1: TEMPL[f'grid{h}x{w}']=[(r,c) for r in range(0,28,h) for c in range(0,34,w)]
for w in range(2,9):
    for k in range(1,w):
        TEMPL[f'diag{w}_{k}']=[(r,c) for r in range(28) for c in range(34) if (c-k*r)%w==0]
def mirror(T):
    m=max(c for r,c in T);return [(r,m-c) for r,c in T]
for k in list(TEMPL): TEMPL[k+'_mirror']=mirror(TEMPL[k])
def tabcoinc(ar,ac,mir):
    S=set()
    for i,(L,c) in enumerate(cells):
        r,cc=L-ar,c-ac
        if 0<=r<len(TAB):
            row=TAB[r][::-1] if mir else TAB[r]
            if 0<=cc<len(row) and row[cc]==CT0[i]: S.add(i)
    return S
def select(T,ar,ac):
    H=set((ar+r,ac+c) for r,c in T);return set(i for i,cl in enumerate(cells) if cl in H)
# ---- selections
SEL={}
for an,(ar,ac) in ANCHORS.items():
    for tn,T in TEMPL.items():
        S=select(T,ar,ac)
        if 0<len(S)<N: SEL[(an,tn)]=S
    for mir in (False,True):
        S=tabcoinc(ar,ac,mir)
        if 0<len(S)<N: SEL[(an,'tableau_coincide',mir)]=S
# text-derived: K1/K2/K3 plaintext aligned to K4 by index (start and end) and plaintext panel overlay at anchor offsets
SETS={'KRYPTOS':set('KRYPTOS'),'YAR':set('YAR'),'QUA':set('QUA'),'LOE':set('LOE'),'X':{'X'},'E':{'E'}}
for ch in ST: SETS.setdefault(ch,{ch})
for tn,T in (('K1P',K1P),('K2P',K2P),('K3P',K3P)):
    for al in ('start','end'):
        s=(T*3)[:N] if al=='start' else (T*3)[-N:]
        for sn,st in SETS.items(): SEL[('text',tn,al,sn)]=set(i for i in range(N) if s[i] in st)
for an,(ar,ac) in ANCHORS.items():
    for sn,st in SETS.items():
        S=set(i for i,(L,c) in enumerate(cells) if 0<=L-ar<len(PB) and 0<=c-ac<len(PB[L-ar]) and PB[L-ar][c-ac] in st)
        if 0<len(S)<N: SEL[(an,'PBoverlay',sn)]=S
print('selections',len(SEL),flush=True)
# ---- permutations: holes first / holes last (raster order), both directions
P=[];L=[]
for k,S in SEL.items():
    a=[i for i in range(N) if i in S];b=[i for i in range(N) if i not in S]
    for nm,o in (('holes_first',a+b),('holes_last',b+a)):
        o=np.array(o);P+=[o,np.argsort(o)];L+=[(k,nm,'a'),(k,nm,'b')]
# Josephus / skip readings with step sequences
STEPS={'morseE':[2,1,6,3,2,2,5,3,1,1],'morseE_noLast':[2,1,6,3,2,2,5,3,1],'coords':[3,8,5,7,6,5,7,7,8,4,4],
       'YAR0':[24,0,17],'YAR1':[25,1,18],'misspell_QUA0':[16,20,0],'misspell_LOE0':[11,14,4],'KRYPTOS':[ST.index(x) for x in 'KRYPTOS']}
for sn,seq in STEPS.items():
    for plus in (0,1):
        for st in range(N):
            ring=list(range(N));pos=st;o=[];j=0
            while ring:
                pos=(pos+seq[j%len(seq)]+plus)%len(ring) if o or plus or True else pos
                o.append(ring.pop(pos));j+=1
                if ring: pos%=len(ring)
            o=np.array(o);P+=[o,np.argsort(o)];L+=[('josephus',sn,plus,st,'a'),('josephus',sn,plus,st,'b')]
P=np.array(P,dtype=np.int16);u,ix=np.unique(P,axis=0,return_index=True);ix=sorted(ix);P=P[ix];L=[L[i] for i in ix]
print('unique perms',len(P),flush=True)
if __name__=='__main__':
    t=time.time();SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for a,ch in ANCH.items(): pt[a]=ch
    i=[x for x,l in enumerate(L) if l[0][:2]==('Y','grid3x5')][0];m=P[i];ct=['?']*N
    for a in range(N): ct[m[a]]=KA[(KA.index(pt[a])+KA.index('MASKYAR'[a%7]))%26]
    print('control',[h for h in run(''.join(ct),ANCH,P[i:i+1],[L[i]],Pmax=12) if h[0]=='key'][:1],flush=True)
    for nm,x in (('K4',CT0),('null',SH)):
        H=[];qq=0
        for c0 in range(0,len(P),20000):
            H+=run(x,ANCH,P[c0:c0+20000],L[c0:c0+20000],Pmax=24);qq+=qii(x,P[c0:c0+20000],24)
        print(nm,'periodic hits',len(H),[h for h in H if h[0]=='key'][:3],'any',sum(1 for h in H if h[0]=='any'),'| QI/II',qq,round(time.time()-t),flush=True)
    # tableau overlay as key at every offset (straight and mirrored), Vig/Beau/varBeau, ST/KA
    best=[]
    for x,nm in ((CT0,'K4'),(SH,'null')):
        bb=(0,)
        for ar in range(-27,28):
            for ac in range(-34,35):
                for mir in (False,True):
                    ks={}
                    for i in CR:
                        Lr,c=cells[i];r,cc=Lr-ar,c-ac
                        if 0<=r<len(TAB):
                            row=TAB[r][::-1] if mir else TAB[r]
                            if 0<=cc<len(row): ks[i]=row[cc]
                    if len(ks)<18: continue
                    for A in (ST,KA):
                        for rn in ('v','b','vb'):
                            h=0
                            for i,k in ks.items():
                                cc_,kk=A.index(x[i]),A.index(k)
                                v=(cc_-kk)%26 if rn=='v' else (kk-cc_)%26 if rn=='b' else (cc_+kk)%26
                                h+=A[v]==ANCH[i]
                            bb=max(bb,(h,len(ks),ar,ac,mir,A[:3],rn))
        print(nm,'tableau-overlay key best',bb,flush=True)
    print('anchors',ANCHORS)
