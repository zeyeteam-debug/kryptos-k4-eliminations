# Re-run every ST/KA-only test family with all 268 sculpture-keyed alphabets (+ST, KA), direct (no transposition).
import numpy as np,itertools,random,time,sys
from stack import CT0,ANCH,solve_mod
from sculptalpha import ALS as SALS
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
ALS=dict(SALS);ALS['ST']=ST;ALS['KA']=KA
AN=np.array(sorted(ANCH))
def Kmat(ct):
    labs=[];rows=[]
    for nm,A in ALS.items():
        C=np.array([A.index(ct[i]) for i in AN]);P=np.array([A.index(ANCH[i]) for i in AN])
        for rule,k in (('vig',(C-P)%26),('beau',(C+P)%26)):
            labs.append((nm,rule));rows.append(k)
    return labs,np.array(rows)
def nullspace_left(A,m):
    from mech import nullspace_left as nl;return nl(A,m)
# ---------- A: linear models ----------
def linear_models():
    import mech
    M=[(l,A) for l,A in mech.ALL]
    ratios={'Metonic':(235,19),'inv':(19,235),'sid':(254,19),'Saros':(223,18),'drac':(242,223),'anom':(239,223),'Callip':(940,76),'exel':(669,54)}
    for nm,(p,q) in ratios.items():
        for s in range(q): M.append((('ratio',nm,s),np.array([AN,(p*AN+s)//q,np.ones(24,int)]).T))
    for s in range(1,24):
        for t in range(24): M.append((('ring',s,t),np.array([(s*AN+t)%24,np.ones(24,int)]).T))
    return M
def linear_test(K,models):
    hits=[]
    for lab,A in models:
        for m in (2,13):
            pass
        n2=nullspace_left(A,2);n13=nullspace_left(A,13);red=len(n13)
        if red<3: continue
        ok=np.ones(len(K),bool)
        if len(n2): ok&=np.all((K%2)@n2.T%2==0,axis=1)
        ok&=np.all((K%13)@n13.T%13==0,axis=1)
        for j in np.nonzero(ok)[0]: hits.append((red,lab,j))
    return hits
# ---------- B: safe piecewise slopes (3 phases, slopes +-1..3), direct ----------
def safe_test(K):
    SL=[-3,-2,-1,1,2,3];ps=AN[None,:];hits=[]
    EQ={s:((K-s*ps)%26)[:,1:]==((K-s*ps)%26)[:,:-1] for s in SL}
    PRE={s:np.argmin(np.c_[EQ[s],np.zeros(len(K),bool)],axis=1) for s in SL}
    SUF={s:np.argmin(np.c_[EQ[s][:,::-1],np.zeros(len(K),bool)],axis=1) for s in SL}
    for s1,s2,s3 in itertools.product(SL,repeat=3):
        x=PRE[s1]+1;y=23-SUF[s3];ok=x>=y
        for j in np.nonzero(~ok)[0]:
            if np.all(EQ[s2][j,x[j]:max(y[j]-1,x[j])]): ok[j]=True
        for j in np.nonzero(ok)[0]: hits.append(((s1,s2,s3),j))
    return hits
# ---------- C: clock-state key bits (Mengenlehreuhr / bars / DCF77) ----------
def mask_test(K):
    from clocks import masks
    hits=[];n=0
    for name,(b,r) in masks.items():
        L=len(b)
        for o in range(L):
            j=(AN+o)%L;lab=np.where(b[j]==1,1+r[j],0);gs=np.unique(lab)
            if 24-len(gs)<4: continue
            ok=np.ones(len(K),bool);n+=1
            for g in gs:
                idx=np.nonzero(lab==g)[0]
                if len(idx)>1: ok&=np.all(K[:,idx]==K[:,idx[:1]],axis=1)
                if not ok.any(): break
            for jj in np.nonzero(ok)[0]: hits.append((name,o,jj))
    return hits,n
# ---------- D: running key, exact substring search ----------
def texts():
    from rk_trans import load_texts
    T,big=load_texts()
    T['carter']=''.join(c for c in open('texts/carter_vol1.txt',errors='ignore').read().upper() if 'A'<=c<='Z')
    from wcity import streams
    for k,v in streams().items(): T['wz'+str(k)]=v+v
    return T
def rk_test(ct,T):
    hits=[];g=63-21
    for nm,A in ALS.items():
        C=[A.index(ct[i]) for i in AN];P=[A.index(ANCH[i]) for i in AN]
        for rule in ('vig','beau','var'):
            ki=[((c-p) if rule=='vig' else (c+p) if rule=='beau' else (p-c))%26 for c,p in zip(C,P)]
            for conv,KAL in (('A',A),('ST',ST)):
                key=''.join(KAL[x] for x in ki);f1=key[:13];f2=key[13:]
                for tn,txt in T.items():
                    for rev in (0,1):
                        s=txt[::-1] if rev else txt;st=s.find(f1)
                        while st!=-1:
                            if s[st+g:st+g+11]==f2: hits.append((nm,rule,conv,tn,rev,st))
                            st=s.find(f1,st+1)
    return hits
# ---------- E: Hill / affine-Hill ----------
def hill_test(ct):
    from hillx import test
    hits=[]
    for nm,A in ALS.items():
        for n in (2,3,4,5):
            for o in range(n):
                for L in range(0,7):
                    r=test(ct,n,o,A,L)
                    if r and r[0] and r[1]>=2: hits.append((nm,n,o,L,r[1]))
    return hits
# ---------- F: autokey ----------
def autokey_test(ct):
    hits=[];an=dict(ANCH)
    for nm,A in ALS.items():
        for L in range(1,41):
            for kind in ('pt','ct'):
                for rule in ('vig','beau','var'):
                    ok=True;chk=0
                    for i in AN:
                        src=an.get(i-L) if kind=='pt' else (ct[i-L] if i-L>=0 else None)
                        if src is None: continue
                        p=A.index(an[i]);k=A.index(src);c=A.index(ct[i])
                        e=[(p+k)%26,(k-p)%26,(p-k)%26][['vig','beau','var'].index(rule)];chk+=1
                        if e!=c: ok=False;break
                    if ok and chk>=4: hits.append((nm,L,kind,rule,chk))
    return hits
if __name__=='__main__':
    t=time.time();SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    print('alphabets',len(ALS),flush=True)
    M=linear_models();print('linear models',len(M),flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        labs,K=Kmat(ct)
        h=linear_test(K,M);print('A linear',nm,len(h),[(r,l,labs[j]) for r,l,j in h[:3]],round(time.time()-t),flush=True)
        h=safe_test(K);print('B safe',nm,len(h),[(s,labs[j]) for s,j in h[:3]],round(time.time()-t),flush=True)
        h,n=mask_test(K);print('C masks',nm,'configs',n,'x',len(K),'hits',len(h),h[:3],round(time.time()-t),flush=True)
        print('E hill',nm,hill_test(ct)[:3],round(time.time()-t),flush=True)
        print('F autokey',nm,autokey_test(ct)[:3],round(time.time()-t),flush=True)
    T=texts();print('texts',len(T),sum(len(v) for v in T.values()),round(time.time()-t),flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        h=rk_test(ct,T);print('D running key',nm,len(h),h[:3],round(time.time()-t),flush=True)
