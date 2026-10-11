# Clock-state masks: Mengenlehreuhr (1440 x seconds lamp), fill-bar clocks (water clock levels), DCF77 minute frames.
import numpy as np,itertools,random,datetime
from stack import CT0,ANCH
from gen import run
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
AN=np.array(sorted(ANCH));SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
def mengen(h,m,s):
    rows=[[s],[1]*(h//5)+[0]*(4-h//5),[1]*(h%5)+[0]*(4-h%5),[1]*(m//5)+[0]*(11-m//5),[1]*(m%5)+[0]*(4-m%5)]
    bits=[b for r in rows for b in r];rowid=[k for k,r in enumerate(rows) for _ in r]
    return np.array(bits),np.array(rowid)
masks={}  # name -> (bits[L], rowid[L])
for h in range(24):
    for m in range(60):
        for s in (0,1):
            b,r=mengen(h,m,s);masks[('mengen',h,m,s)]=(b,r);masks[('mengenR',h,m,s)]=(b[::-1],r[::-1])
for L in (12,24,60):
    for f in range(1,L): masks[('bar',L,f)]=(np.array([1]*f+[0]*(L-f)),np.zeros(L,int))
def dcf(dt):
    b=[0]*59
    b[17],b[18]=(1,0) if 3<dt.month<10 else (0,1)  # crude CEST/CET
    b[20]=1
    def bcd(v,pos,n):
        d=[1,2,4,8,10,20,40,80]
        for i in range(n):
            pass
        bits=[];t=v%10;u=v//10
        for k in range(4): bits.append((t>>k)&1)
        for k in range(n-4): bits.append((u>>k)&1)
        for i,x in enumerate(bits): b[pos+i]=x
        return bits
    p=bcd(dt.minute,21,7);b[28]=sum(p)%2
    p=bcd(dt.hour,29,6);b[35]=sum(p)%2
    q=bcd(dt.day,36,6)+[ (dt.isoweekday()>>k)&1 for k in range(3)]
    for i in range(3): b[42+i]=q[6+i]
    q+=bcd(dt.month,45,5);q+=bcd(dt.year%100,50,8);b[58]=sum(q)%2
    return np.array(b)
for d in (datetime.date(1990,11,3),datetime.date(1989,11,9),datetime.date(1990,10,3),datetime.date(1988,11,3)):
    for h in range(24):
        for m in range(60):
            dt=datetime.datetime(d.year,d.month,d.day,h,m)
            masks[('dcf77',str(d),h,m)]=(dcf(dt),np.zeros(59,int))
print('masks',len(masks),flush=True)
def Kmat(ct):
    out={}
    for an,A in (('ST',ST),('KA',KA)):
        C=np.array([A.index(ct[i]) for i in AN]);P=np.array([A.index(ANCH[i]) for i in AN])
        out[(an,'vig')]=(C-P)%26;out[(an,'beau')]=(C+P)%26
    return out
# (1) lamp states as key: k = c + sum_row w_row * lamp   (unlit all share c; lit lamps share value per row)
def keytest(ct,rowweights=True):
    Ks=Kmat(ct);hits=[];n=0
    for name,(b,r) in masks.items():
        L=len(b)
        for o in range(L):
            j=(AN+o)%L;lab=np.where(b[j]==1,1+(r[j] if rowweights else 0),0)
            ng=len(set(lab.tolist()));red=24-ng
            if red<4: continue
            for vk,K in Ks.items():
                n+=1;ok=True
                for g in set(lab.tolist()):
                    if len(set(K[lab==g].tolist()))>1: ok=False;break
                if ok: hits.append((red,name,o,vk))
    return hits,n
# (2) grille transposition: in each L-block read lit cells first, then unlit (and reverse); also global (whole text)
def grille_perms():
    P=[];Lb=[]
    for name,(b,r) in masks.items():
        L=len(b);tile=np.array([b[i%L] for i in range(97)])
        for mode in ('block','global'):
            if mode=='global': order=list(np.nonzero(tile==1)[0])+list(np.nonzero(tile==0)[0])
            else:
                order=[]
                for s in range(0,97,L):
                    blk=np.arange(s,min(97,s+L));order+=list(blk[tile[blk]==1])+list(blk[tile[blk]==0])
            order=np.array(order)
            P.append(order);Lb.append((name,mode,'a'))
            inv=np.argsort(order);P.append(inv);Lb.append((name,mode,'b'))
    P=np.array(P,dtype=np.int16);u,ix=np.unique(P,axis=0,return_index=True)
    return P[sorted(ix)],[Lb[i] for i in sorted(ix)]
# (3) null masks: lit (or unlit) cells, tiled, are nulls; cribs must be untouched; test reduced text, any periodic key <=24
def nulltest(ct):
    hits=[];n=0;seen=set()
    for name,(b,r) in masks.items():
        L=len(b)
        for inv in (0,1):
            for o in range(L):
                tile=np.array([b[(i+o)%L]^inv for i in range(97)]);nulls=set(np.nonzero(tile==1)[0].tolist())
                if not nulls or len(nulls)>40 or nulls&set(ANCH): continue
                key=tuple(sorted(nulls))
                if key in seen: continue
                seen.add(key);keep=[i for i in range(97) if i not in nulls];ct2=''.join(ct[i] for i in keep)
                an2={keep.index(i):c for i,c in ANCH.items()};n+=1
                H=run(ct2,an2,np.arange(len(ct2),dtype=np.int16)[None,:],[(name,inv,o)],Pmax=24)
                if H: hits.append(H)
    return hits,n
if __name__=='__main__':
    # controls
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    b,r=masks[('mengen',14,37,1)]
    w=[3,0,7,0,19,5]
    ctl=''.join(KA[(KA.index(c)+4+(w[r[(i+5)%24]+1] if b[(i+5)%24] else 0))%26] for i,c in enumerate(pt))
    h,_=keytest(ctl);print('control key',[x for x in h if x[1]==('mengen',14,37,1) and x[2]==5][:1],flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        h,n=keytest(ct);print('(1) key',nm,'tested',n,'hits',len(h),sorted(h,key=lambda x:-x[0])[:5],flush=True)
    P,Lb=grille_perms();print('grille perms',len(P),flush=True)
    import mech
    from mech import check,special
    keep_r={12,24,60}|{x for t in special.values() for x in t if x>=2}
    sel=[(l,A) for (l,A) in mech.ALL if (l[0]=='clock' and l[1] in keep_r) or l[0] in('clock3','gear')]
    mech.build(sel);s2=[m for m,r in zip(sel,mech.RED) if r>=6];mech.build(s2)
    for nm,ct in (('K4',CT0),('null',SH)):
        H=[];[H.extend(run(ct,ANCH,P[s:s+20000],Lb[s:s+20000],Pmax=24)) for s in range(0,len(P),20000)]
        h2,e=check(ct,P,Lb,minred=6)
        print('(2) grille',nm,'periodic/any hits',len(H),H[:3],'clock/gear',len(h2),flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        h,n=nulltest(ct);print('(3) null masks',nm,'tested',n,'hits',len(h),h[:3],flush=True)
