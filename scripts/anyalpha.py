"""Alphabet-free test: periodic substitution with P completely arbitrary alphabets (covers any keyword/mixed alphabet,
Vigenere/Beaufort/Quagmire I-IV with any alphabet), optionally after a transposition.
Within one key class: same PT letter <-> same CT letter."""
import numpy as np,itertools,random,sys
from stack import CT0,ANCH
N=97;AN=np.array(sorted(ANCH));PR=np.array(list(itertools.combinations(range(len(AN)),2)))
Pl=np.array([ord(ANCH[i]) for i in AN]);Peq=Pl[PR[:,0]]==Pl[PR[:,1]]
def run(ct,Ms,labs,Pmax=12,strong=4):
    ci=np.array([ord(x) for x in ct]);C=ci[Ms[:,AN]];Ceq=C[:,PR[:,0]]==C[:,PR[:,1]]
    bad=Peq^Ceq;hits=[]
    for P in range(1,Pmax+1):
        # pt-indexed classes
        same=(AN[PR[:,0]]%P)==(AN[PR[:,1]]%P)
        st=int((same&Peq).sum())
        if st>=strong:
            ok=~np.any(bad[:,same],axis=1)
            for x in np.nonzero(ok)[0]: hits.append(('pt',P,st,labs[x]))
        # ct-indexed classes
        R=Ms[:,AN]%P;sm=R[:,PR[:,0]]==R[:,PR[:,1]]
        ok=~np.any(bad&sm,axis=1)&((sm&Peq).sum(axis=1)>=strong)
        for x in np.nonzero(ok)[0]: hits.append(('ct',P,int((sm[x]&Peq).sum()),labs[x]))
    return hits
def families():
    from tauto import perms
    import cyl,counter
    for W in range(2,10):
        M,Mi=perms(W);yield np.concatenate([M,Mi]),[('col',W,'enc',i) for i in range(len(M))]+[('col',W,'dec',i) for i in range(len(M))]
    Ms,l=cyl.family();yield Ms,[('cyl',)+x for x in l]
    Ms,l=counter.family();yield Ms,[('rot',)+x for x in l]
def runall(ct,label):
    tot=[];ident=np.arange(N,dtype=np.int16)[None,:]
    tot+=run(ct,ident,[('none',)])
    for Ms,labs in families():
        for s in range(0,len(Ms),50000): tot+=run(ct,Ms[s:s+50000],labs[s:s+50000])
    return tot
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    from tauto import colmap
    idx=colmap(7,(2,5,0,6,3,1,4));random.seed(11)
    alphs=[random.sample('ABCDEFGHIJKLMNOPQRSTUVWXYZ',26) for _ in range(3)]
    x=''.join(alphs[i%3][ord(c)-65] for i,c in enumerate(pt))   # 3 random alphabets
    ctc=''.join(x[i] for i in idx)
    print('direct, no transposition, which P are possible for K4:',[h[1:3] for h in run(CT0,np.arange(N,dtype=np.int16)[None,:],[('none',)],Pmax=24,strong=1)][:30])
    h=runall(ctc,'ctrl');print('control (col W7 + 3 random alphabets): true found',any(t[3][:2]==('col',7) and t[1]==3 for t in h),'total hits',len(h),flush=True)
    random.seed(5);s=list(CT0);random.shuffle(s);hn=runall(''.join(s),'null');print('null hits',len(hn),flush=True)
    hk=runall(CT0,'K4');print('K4 hits',len(hk));
    from collections import Counter;print(Counter((h[0],h[1]) for h in hk).most_common(10));print(hk[:10])
