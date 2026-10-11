# Fold the sheet in two / four (or stack the panels) so that K1-K3 (solved or unsolved, with or without the
# misspellings, with or without the raised YAR) and the tableau land on K4. The letter landing on each K4 cell is used
# (a) directly as key, (b) as running key + periodic key, (c) as a mask (coincidences) giving reading orders.
exec(open('handmask.py').read().split('# ---- templates')[0].replace("print(","(lambda *a,**k:None)("))
import random,itertools
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
def board_from_stream(stream):
    it=iter(stream);B=[]
    for L,line in enumerate(LINES):
        row=[]
        for col,ch in enumerate(line):
            if (L<24 or (L==24 and col<26)) and ch!='?': row.append(next(it,'?'))
            else: row.append(ch)
        B.append(''.join(row))
    return B
K1c=K1P.replace('IQLUSION','ILLUSION');K2c=K2P.replace('UNDERGRUUND','UNDERGROUND')
K3s=K3P.replace('DESPARATLY','DESPERATLY');K3f=K3P.replace('DESPARATLY','DESPERATELY')
SRC={}
SRC['solved_asis']=PB
SRC['solved_fixed_sub']=board_from_stream(K1c+K2c+K3s)
SRC['solved_fixed_full']=board_from_stream(K1c+K2c+K3f)
# unsolved with the misspelled letters re-encrypted correctly (K1: KA Vigenere PALIMPSEST, K2: ABSCISSA; K3: letter swap in place)
def vig(p,k,i): return KA[(KA.index(p)+KA.index(k[i%len(k)]))%26]
ctfix=[list(l) for l in LINES]
i1=K1P.index('IQLUSION')+1;L,c=LC[i1];ctfix[L][c]=vig('L','PALIMPSEST',i1)
i2=K2P.index('UNDERGRUUND')+7;L,c=LC[len(K1P)+i2];ctfix[L][c]=vig('O','ABSCISSA',i2)
L,c=MIS['DESPARATLY_A'];ctfix[L][c]='E'
assert ''.join(LINES[LC[i1][0]])[LC[i1][1]]=='K'
SRC['unsolved_asis']=LINES;SRC['unsolved_fixed']=[''.join(r) for r in ctfix]
TABS=TAB
def withyar(B,mode):
    if mode=='keep': return B
    B=[list(r) for r in B]
    for c in (3,4,6): B[14][c]='?'
    return [''.join(r) for r in B]
def get(B,r,c):
    if 0<=r<len(B) and 0<=c<len(B[r]):
        ch=B[r][c];return ch if ch.isalpha() else None
    return None
FOLDS={}
for wref in ('sheet33','rowlen'):
    def W(r,B,wref=wref): return 33 if wref=='sheet33' else len(B[r])
    FOLDS[f'half_horizontal']=lambda r,c,B:(27-r,c)
    FOLDS[f'half_vertical_{wref}']=lambda r,c,B,W=W:(r,W(r,B)-1-c)
    FOLDS[f'quarter_{wref}']=lambda r,c,B,W=W:(27-r,W(27-r,B)-1-c)
for dc in range(-8,9):
    FOLDS[f'stack_panels_dc{dc}']=lambda r,c,B,dc=dc:(r-14,c+dc)
    FOLDS[f'stack_panels_flipped_dc{dc}']=lambda r,c,B,dc=dc:(27-r,c+dc)  # = horizontal fold with shift
FOLDS['quarter_into_K4corner']=lambda r,c,B:(r-24 if r>=24 else r, c)  # top-left quarter laid on K4 rows
combos=[]
for sn,B0 in SRC.items():
    for ym in ('keep','blank'):
        B=withyar(B0,ym)
        for fn,f in FOLDS.items():
            combos.append(((sn,ym,fn),[get(B,*f(r,c,B)) for (r,c) in cells]))
for fn in ('straight','mirror31','mirror_rowlen'):
    for dr in range(-4,5):
        ks=[]
        for (r,c) in cells:
            rr=r+dr
            if not 0<=rr<28: ks.append(None);continue
            row=TABS[rr];cc=c if fn=='straight' else (31-c if fn=='mirror31' else len(row)-1-c)
            ks.append(row[cc] if 0<=cc<len(row) else None)
        combos.append((('tableau',fn,dr),ks))
print('combos',len(combos),flush=True)
AN=sorted(ANCH)
def direct(ct,ks):
    best=(0,0,'')
    for A,an in ((ST,'ST'),(KA,'KA')):
        for m in ('v','b','vb'):
            h=n=0
            for i in AN:
                k=ks[i]
                if not k: continue
                n+=1;cc,kk=A.index(ct[i]),A.index(k)
                v=(cc-kk)%26 if m=='v' else (kk-cc)%26 if m=='b' else (cc+kk)%26
                h+=A[v]==ANCH[i]
            best=max(best,(h,n,an+m))
    return best
def withperiodic(ct,ks):
    hits=[]
    for A,an in ((ST,'ST'),(KA,'KA')):
        for m in ('v','b','vb'):
            res={}
            for i in AN:
                k=ks[i]
                if not k: continue
                cc,kk,pp=A.index(ct[i]),A.index(k),A.index(ANCH[i])
                # residual periodic key letter
                res[i]=(cc-kk-pp)%26 if m=='v' else (kk-cc-pp)%26 if m=='b' else (cc+kk-pp)%26
            if len(res)<16: continue
            for p in range(1,27):
                cl={};ok=True
                for i,v in res.items():
                    if cl.setdefault(i%p,v)!=v: ok=False;break
                if ok and len(res)-len({i%p for i in res})>=6: hits.append((an+m,p))
    return hits
for nm,x in (('K4',CT0),('null',SH)):
    D=sorted(((direct(x,ks),lab) for lab,ks in combos),reverse=True)
    Hp=[(lab,h) for lab,ks in combos for h in withperiodic(x,ks)]
    print(nm,'direct best',D[:3],'| +periodic hits',len(Hp),Hp[:3],flush=True)
# control: K4-like text whose key is the solved panel folded horizontally (KA Vigenere) + period 5
lab,ks=[cb for cb in combos if cb[0]==('solved_asis','keep','half_horizontal')][0]
pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
for a,ch in ANCH.items(): pt[a]=ch
ctl=''.join(KA[(KA.index(pt[i])+KA.index(ks[i] or 'A')+KA.index('FOLDS'[i%5]))%26] for i in range(97))
print('control +periodic',[h for h in withperiodic(ctl,ks) if h[1]==5][:1],'direct-only control',direct(''.join(KA[(KA.index(pt[i])+KA.index(ks[i] or 'A'))%26] for i in range(97)),ks))
# (c) masks: holes where the landing letter equals the K4 letter / is blank -> reading orders
P=[];Lb=[]
for lab,ks in combos:
    for rel in ('eq','blank'):
        S=[i for i in range(97) if (ks[i]==CT0[i] if rel=='eq' else ks[i] is None)]
        if 0<len(S)<97:
            a=S;b=[i for i in range(97) if i not in S]
            for o,on in ((a+b,'first'),(b+a,'last')):
                o=np.array(o);P+=[o,np.argsort(o)];Lb+=[(lab,rel,on,'a'),(lab,rel,on,'b')]
P=np.array(P,dtype=np.int16);u,ix=np.unique(P,axis=0,return_index=True);ix=sorted(ix);P=P[ix];Lb=[Lb[i] for i in ix]
for nm,x in (('K4',CT0),('null',SH)):
    print(nm,'mask perms',len(P),'periodic',len(run(x,ANCH,P,Lb,Pmax=24)),'QI/II',qii(x,P,24),flush=True)
