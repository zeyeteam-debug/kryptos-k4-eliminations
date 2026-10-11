# Length-masking hypothesis: letters inserted/deleted between the cribs (or before ENE) so the key index of the
# BERLINCLOCK block is shifted by D relative to EASTNORTHEAST. Exact tests 1-4 with planted controls + shuffled nulls.
import numpy as np, itertools, random, time, sys
from stack import CT0, ANCH
from gen import colperms
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
AN=np.array(sorted(ANCH));B2=(AN>=63).astype(int);PT=np.array([ord(ANCH[i]) for i in AN])
DS=[d for d in range(-10,11) if d!=0];T0=time.time()
def log(*a): print(*a,round(time.time()-T0),flush=True)
PAIRS=np.array(list(itertools.combinations(range(24),2)));PEQ=PT[PAIRS[:,0]]==PT[PAIRS[:,1]]
def eidx(D): return AN+D*B2
def classes(D,P):
    R=eidx(D)%P;a=[];b=[]
    for r in np.unique(R):
        g=np.nonzero(R==r)[0]
        for x in g[1:]: a.append(g[0]);b.append(x)
    sm=R[PAIRS[:,0]]==R[PAIRS[:,1]]
    return np.array(a,int),np.array(b,int),sm
def Kmats(CL):
    """CL: n x 24 char codes of ct at crib positions -> dict (alph,rule)->n x 24 key values (var == -vig, same consistency)"""
    out={}
    for an,A in (('ST',ST),('KA',KA)):
        lut=np.zeros(128,np.int16);lut[[ord(c) for c in A]]=np.arange(26)
        C=lut[CL];P=lut[PT]
        out[(an,'vig')]=((C-P)%26).astype(np.int8);out[(an,'beau')]=((C+P)%26).astype(np.int8)
    return out
def periodic(CL,Pmax,Ds=DS,anyalpha=True,keep=None):
    """returns list of hits (kind, D, P, alph/rule, redundancy-or-strong, row)"""
    KM=Kmats(CL);hits=[]
    CEQ=(CL[:,PAIRS[:,0]]==CL[:,PAIRS[:,1]]) if anyalpha else None
    for D in Ds:
        for P in range(1,Pmax+1):
            a,b,sm=classes(D,P);red=len(a)
            if red>0:
                for k,K in KM.items():
                    ok=np.all(K[:,b]==K[:,a],axis=1)
                    for x in np.nonzero(ok)[0]: hits.append(('key',D,P,k,red,int(x)))
            if anyalpha:
                strong=int((sm&PEQ).sum())
                if strong>=1:
                    cols=np.nonzero(sm)[0]
                    ok=np.all(CEQ[:,cols]==PEQ[cols],axis=1)
                    for x in np.nonzero(ok)[0]: hits.append(('any',D,P,'-',strong,int(x)))
    return hits
def shuffled(n,seed):
    r=random.Random(seed);out=[]
    for _ in range(n):
        s=list(CT0);r.shuffle(s);out.append(''.join(s))
    return out
def CLof(ct): return np.array([[ord(ct[i]) for i in AN]],np.int16)
# ---------------- controls helpers ----------------
def plant_pt(seed):
    r=random.Random(seed);pt=[r.choice(ST) for _ in range(97)]
    for i,c in ANCH.items(): pt[i]=c
    return pt
def enc_periodic(pt,D,key,A,rule,pos_of=None):
    # key index of pt position i: i (+D if i>=63); pos_of maps pt index -> ct index (transposition)
    ct=['?']*97
    for i,p in enumerate(pt):
        e=i+(D if i>=63 else 0);k=A.index(key[e%len(key)]);pp=A.index(p)
        c=(pp+k)%26 if rule=='vig' else (k-pp)%26
        ct[i if pos_of is None else pos_of[i]]=A[c]
    return ''.join(ct)
def enc_anyalpha(pt,D,P,seed,pos_of=None):
    r=random.Random(seed);alph=[]
    for _ in range(P): s=list(ST);r.shuffle(s);alph.append(''.join(s))
    ct=['?']*97
    for i,p in enumerate(pt):
        e=i+(D if i>=63 else 0);ct[i if pos_of is None else pos_of[i]]=alph[e%P][ST.index(p)]
    return ''.join(ct)
def summarize(h,thr_key=3,thr_any=4):
    k=[x for x in h if x[0]=='key' and x[4]>=thr_key];a=[x for x in h if x[0]=='any' and x[4]>=thr_any]
    return k,a
# ================= TEST 1: direct periodic =================
def test1(NS=2000):
    log('=== T1 direct periodic, D=-10..10 (D!=0), P=1..26, ST/KA vig/beau(var==vig), + arbitrary alphabets')
    # controls
    pt=plant_pt(1);ct=enc_periodic(pt,4,'PALIMPS',KA,'vig');h=periodic(CLof(ct),26)
    log('ctrl key KA vig P7 D+4 recovered:',any(x[1]==4 and x[2]==7 and x[3]==('KA','vig') for x in h),
        'all consistent red>=6:',sorted({(x[1],x[2],x[3]) for x in h if x[0]=='key' and x[4]>=6})[:6])
    ct=enc_anyalpha(pt,-6,5,3);h=periodic(CLof(ct),26)
    log('ctrl anyalpha P5 D-6 recovered:',any(x[0]=='any' and x[1]==-6 and x[2]==5 for x in h),
        'strong of it:',[x[4] for x in h if x[0]=='any' and x[1]==-6 and x[2]==5])
    hk=periodic(CLof(CT0),26)
    shuf=shuffled(NS,11);CLn=np.array([[ord(s[i]) for i in AN] for s in shuf],np.int16)
    hn=periodic(CLn,26)
    for thr in (2,3,4,6):
        kk=[x for x in hk if x[0]=='key' and x[4]>=thr];kn=[x for x in hn if x[0]=='key' and x[4]>=thr]
        log(f' key red>={thr}: K4 {len(kk)} | null mean/shuffle {len(kn)/NS:.3f} | null frac shuffles with >=K4: {np.mean(np.bincount([x[5] for x in kn],minlength=NS)>=len(kk)):.3f}')
    log(' K4 key hits red>=3 (D,P,alph,rule,red; D%P==0 => baseline-equivalent):',
        sorted([(x[1],x[2],x[3],x[4],'BASE' if x[1]%x[2]==0 else '') for x in hk if x[0]=='key' and x[4]>=3],key=lambda t:-t[3])[:25])
    for thr in (2,4,6,8):
        ka=[x for x in hk if x[0]=='any' and x[4]>=thr];na=[x for x in hn if x[0]=='any' and x[4]>=thr]
        log(f' anyalpha strong>={thr}: K4 {len(ka)} | null mean {len(na)/NS:.3f} | frac null >=K4: {np.mean(np.bincount([x[5] for x in na],minlength=NS)>=len(ka)):.3f}')
    log(' K4 anyalpha strong>=4 non-baseline:',sorted([(x[1],x[2],x[4]) for x in hk if x[0]=='any' and x[4]>=4 and x[1]%x[2]],key=lambda t:-t[2])[:20])
# ================= TEST 2: transpositions, pt-indexed key =================
def perms():
    from k3rot import compose
    M=[];labs=[]
    for W in range(2,9):
        c=colperms(97,W);M.append(c.astype(np.int16));n=len(c)//2
        labs+=[('col',W,j,'w') for j in range(n)]+[('col',W,j,'inv') for j in range(n)]
    R=[];rl=[]
    for w1 in range(2,97):
        for w2 in range(2,97):
            for d1 in (True,False):
                for d2 in (True,False):
                    idx=compose(97,w1,d1,w2,d2);R.append(np.argsort(idx));rl.append(('k3',w1,d1,w2,d2,'fwd'))
                    R.append(idx);rl.append(('k3',w1,d1,w2,d2,'inv'))
    M.append(np.array(R,np.int16));labs+=rl
    return np.concatenate(M),labs
def test2():
    log('=== T2 transpositions + pt-indexed periodic key (P<=12) with block-2 key index +D')
    Ms,labs=perms();log('perms',len(Ms))
    # control: plant via a K3 double rotation and a columnar
    pt=plant_pt(5)
    for j in (labs.index(('k3',24,True,8,False,'fwd')),labs.index(('col',7,1234,'inv'))):
        ct=enc_periodic(pt,-3,'KRYPTOSAB'[:9],ST,'beau',pos_of=Ms[j]);CL=np.array([[ord(ct[Ms[j][i]]) for i in AN]],np.int16)
        h=periodic(CL,12)
        log(' ctrl key',labs[j],'ST beau P9 D-3 recovered:',any(x[1]==-3 and x[2]==9 and x[3]==('ST','beau') for x in h))
        ct=enc_anyalpha(pt,7,4,9,pos_of=Ms[j]);CL=np.array([[ord(ct[Ms[j][i]]) for i in AN]],np.int16)
        h=periodic(CL,12);log(' ctrl anyalpha',labs[j],'P4 D+7 recovered:',any(x[0]=='any' and x[1]==7 and x[2]==4 for x in h))
    res={}
    for nm,ct in (('K4',CT0),('null1',shuffled(1,21)[0]),('null2',shuffled(1,22)[0])):
        cl=np.array([ord(c) for c in ct],np.int16);CL=cl[Ms[:,AN]]
        h=periodic(CL,12);res[nm]=h
        k=[x for x in h if x[0]=='key'];a=[x for x in h if x[0]=='any' and x[4]>=4]
        kb=[x for x in k if x[1]%x[2]];ab=[x for x in a if x[1]%x[2]]
        log(f' {nm}: key-consistent {len(k)} (non-baseline {len(kb)}), min red {min([x[4] for x in k],default=None)}; anyalpha strong>=4 {len(a)} (non-baseline {len(ab)})')
        for x in sorted(kb,key=lambda t:-t[4])[:5]: log('   key',x[:5],labs[x[5]])
        for x in sorted(ab,key=lambda t:-t[4])[:5]: log('   any',x[:5],labs[x[5]])
# ================= TEST 3: running key with gap shift =================
def test3():
    from allalpha import texts
    from sculptalpha import ALS as SA
    ALS=dict(SA);ALS['ST']=ST;ALS['KA']=KA
    T=texts();big='|'.join(T.values());rev=big[::-1];log('=== T3 running key; texts',len(T),'letters',len(big),'alphabets',len(ALS))
    def keys(ct):
        out=[]
        for nm,A in ALS.items():
            C=[A.index(ct[i]) for i in AN];P=[A.index(ANCH[i]) for i in AN]
            for rule in ('vig','beau','var'):
                ki=[((c-p) if rule=='vig' else (c+p) if rule=='beau' else (p-c))%26 for c,p in zip(C,P)]
                for conv,KAL in (('A',A),('ST',ST)):
                    out.append(((nm,rule,conv),''.join(KAL[x] for x in ki)))
        return out
    def search(ct,ks=None):
        hits=[];n1=n2=0
        for lab,key in (ks or keys(ct)):
            f1=key[:13];f2=key[13:]
            for d,s in (('fwd',big),('rev',rev)):
                st=s.find(f1)
                while st!=-1:
                    n1+=1
                    for D in DS:
                        if s[st+42+D:st+53+D]==f2: hits.append((lab,d,D,st))
                    st=s.find(f1,st+1)
                st=s.find(f2)
                while st!=-1: n2+=1;st=s.find(f2,st+1)
        return hits,n1,n2
    # control: key = Carter text, block2 key taken at offset +42+D with D=+5; alphabet from ALS
    A=ALS['PALIMPSEST/ST'];o=123456;src=big;pt=plant_pt(8);D=5
    ks=src[o:o+13]+src[o+42+D:o+53+D];ct=list('?'*97)
    for j,i in enumerate(AN): ct[i]=A[(A.index(ANCH[i])+A.index(ks[j]))%26]
    ct=''.join(c if c!='?' else 'A' for c in ct)
    h,_,_=search(ct,[x for x in keys(ct) if x[0][0]=='PALIMPSEST/ST'])
    log(' ctrl D+5 PALIMPSEST/ST vig recovered:',[x for x in h if x[2]==5][:2])
    for nm,ct in (('K4',CT0),('null',shuffled(1,31)[0])):
        h,n1,n2=search(ct);log(f' {nm}: gap-shift hits {len(h)} {h[:5]} | 13-letter ENE-key fragment occurrences {n1}, 11-letter BC-key fragment occurrences {n2}')
# ================= TEST 4: linear clock/gear models with shifted block-2 positions =================
def models_for(E):
    def gear(ps): return np.array([(E%p==r).astype(int) for p in ps for r in range(p)]).T
    def clock(r,s): return np.array([E,(E+s)//r,np.ones(24,int)]).T
    def clock3(r1,r2,s): return np.array([E,(E+s)//r1,(E+s)//(r1*r2),np.ones(24,int)]).T
    import mech
    M=[]
    for p1 in range(2,31):
        for p2 in range(p1+1,31): M.append((('gear',p1,p2),gear((p1,p2))))
    for nm,t in mech.special.items():
        ps=tuple(sorted(set(x if x>0 else 26 for x in t)-{1}))
        if ps: M.append((('gear',nm)+ps,gear(ps)))
    for r in range(2,61):
        for s in range(r): M.append((('clock',r,s),clock(r,s)))
    for r1,r2 in ((12,12),(24,12),(12,5),(60,12),(5,12)):
        for s in range(r1): M.append((('clock3',r1,r2,s),clock3(r1,r2,s)))
    return M
def test4(NS=300):
    from mech import nullspace_left as NL
    import mech
    assert len(models_for(AN))==len(mech.ALL)
    log('=== T4 linear models',len(mech.ALL),'x',len(DS),'shifts; direct ST/KA vig/beau; red13>=4')
    def Krows(cts):
        out=[];lab=[]
        for t,ct in enumerate(cts):
            for k,K in Kmats(CLof(ct)).items(): out.append(K[0].astype(np.int64));lab.append((t,k))
        return np.array(out),lab
    pt=plant_pt(4);r,s,D=17,5,-7;key=[0]*97
    for i in range(97):
        e=i+(D if i>=63 else 0);key[i]=(3*e+11*((e+s)//r)+4)%26
    ctc=''.join(KA[(KA.index(p)+key[i])%26] for i,p in enumerate(pt))
    shuf=shuffled(NS,41);Kall,lab=Krows([CT0,ctc]+shuf)
    hits=[];exp=0.0;cache={}
    for D in DS:
        for ml,A in models_for(eidx(D)):
            n2=NL(A,2);n13=NL(A,13);red=len(n13)
            if red<4: continue
            ok=np.ones(len(Kall),bool)
            if len(n2): ok&=np.all((Kall%2)@n2.T%2==0,axis=1)
            ok&=np.all((Kall%13)@n13.T%13==0,axis=1)
            exp+=4*2.0**-len(n2)*13.0**-red
            for x in np.nonzero(ok)[0]: hits.append((lab[x][0],D,ml,lab[x][1],red))
    log(' ctrl clock(17,5) D-7 KA vig recovered:',any(h[0]==1 and h[1]==-7 and h[2]==('clock',17,5) and h[3]==('KA','vig') for h in hits))
    k4=[h for h in hits if h[0]==0];nl=[h for h in hits if h[0]>=2]
    log(f' K4 hits {len(k4)} {k4[:8]} | null {len(nl)} over {NS} shuffles = {len(nl)/NS:.3f}/text | expected/text {exp:.4f}')
if __name__=='__main__':
    which=sys.argv[1:] or ['1','2','3','4']
    for w in which: {'1':test1,'2':test2,'3':test3,'4':test4}[w]()
    log('done')
