# Period 19 (and 38) with keyword-mixed alphabets: Quagmire I/II/III over 370k dictionary words (+ sculpture words),
# Quagmire IV over all pairs of sculpture-derived alphabets. Direct (no transposition).
import numpy as np,time,random
from stack import CT0,ANCH
from sculptalpha import ALS as SALS
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
AN=np.array(sorted(ANCH))
def pairs(P):
    return [(a,b) for ai,a in enumerate(AN) for b in AN[ai+1:] if (b-a)%P==0]
PR={P:np.array(pairs(P)) for P in (19,38)}
def posmat(words,base):
    M=np.zeros((len(words),26),np.int8)
    for r,w in enumerate(words):
        seen=[];
        for c in w+base:
            if c not in seen: seen.append(c)
        for i,c in enumerate(seen): M[r,ord(c)-65]=i
    return M
def check(ct,Mp,Mc,P):
    # Mp: pos of letter in pt alphabet (n x 26), Mc: in ct alphabet. shift s_i = Mc[c_i]-Mp[p_i] (vig) or Mc[c_i]+Mp[p_i] (beau)
    c=np.array([ord(ct[i])-65 for i in AN]);p=np.array([ord(ANCH[i])-65 for i in AN])
    idx={a:j for j,a in enumerate(AN)};pr=PR[P];a=np.array([idx[x] for x in pr[:,0]]);b=np.array([idx[x] for x in pr[:,1]])
    out={}
    for rule in ('vig','beau'):
        S=(Mc[:,c].astype(int)-Mp[:,p]) if rule=='vig' else (Mc[:,c].astype(int)+Mp[:,p]);S%=26
        out[rule]=np.nonzero(np.all(S[:,a]==S[:,b],axis=1))[0]
    return out,len(pr)
if __name__=='__main__':
    t=time.time();SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    import dictsearch as D
    W=sorted(set(D.WORDS)|{a.split('/')[0] for a in SALS});print('words',len(W),'pairs per period',{P:len(v) for P,v in PR.items()},flush=True)
    bases={'ST':ST,'STrev':ST[::-1],'KA':KA}
    # control: Quagmire III, alphabet keyed CLEPSYDRA, period-19 key
    import random as R;rr=R.Random(1);pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,ch in ANCH.items(): pt[i]=ch
    A=''.join(dict.fromkeys('CLEPSYDRA'+ST));key=[rr.randrange(26) for _ in range(19)]
    ctl=''.join(A[(A.index(ch)+key[i%19])%26] for i,ch in enumerate(pt))
    M=posmat(['CLEPSYDRA','KRYPTOS'],ST);o,_=check(ctl,M,M,19);print('control Q3 CLEPSYDRA p19:',{k:list(v) for k,v in o.items()},flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        tot={}
        for bn,base in bases.items():
            M=posmat(W,base);STm=posmat(['A'],ST)[[0]*len(W)];KAm=posmat(['K'],KA)[[0]*len(W)]
            for P in (19,38):
                for qn,(Mp,Mc) in {'Q3':(M,M),'Q1':(M,STm),'Q2':(STm,M),'Q1ka':(M,KAm),'Q2ka':(KAm,M)}.items():
                    o,npairs=check(ct,Mp,Mc,P)
                    for rule,h in o.items():
                        if len(h): tot[(bn,P,qn,rule)]=[W[i] for i in h[:5]]
        print(nm,'Q1/Q2/Q3 dictionary, periods 19/38: consistent ->',tot if tot else 0,round(time.time()-t),flush=True)
        # Quagmire IV over sculpture alphabets x sculpture alphabets
        names=list(SALS);Ms=np.array([[SALS[n].index(chr(65+i)) for i in range(26)] for n in names],np.int8)
        hits=[]
        for j in range(len(names)):
            Mc=np.repeat(Ms[[j]],len(names),axis=0)
            for P in (19,38):
                o,_=check(ct,Ms,Mc,P)
                for rule,h in o.items(): hits+=[(names[i],names[j],P,rule) for i in h]
        print(nm,'Q4 sculpture x sculpture alphabets:',len(hits),hits[:5],round(time.time()-t),flush=True)
