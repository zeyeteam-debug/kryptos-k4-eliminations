import numpy as np,random,collections
from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
TARGETS=['3857','385765','7708','770844','1990','1989','1986','0309','1103','0911','1109','9111989','2111990','3111990','1912','1922','67','675','6750']
def digits(ct,A,base): return ''.join(str(A.index(c)+base) for c in ct)
def feats(ct):
    out={}
    for an,A in (('ST',ST),('KA',KA)):
        for b in (0,1):
            d=digits(ct,A,b);v=np.array([A.index(c)+b for c in ct])
            out[(an,b,'targets')]=[t for t in TARGETS if t in d]
            out[(an,b,'sum')]=int(v.sum());out[(an,b,'sum%26')]=int(v.sum()%26)
            diff=np.diff(v)%26;out[(an,b,'maxrun_equal_diff')]=max(len(list(g)) for _,g in __import__('itertools').groupby(diff))
            out[(an,b,'mod10_chi')]=round(float(((np.bincount(v%10,minlength=10)-9.7)**2/9.7).sum()),1)
            r31=[int(v[i:i+31].sum()) for i in (4,35,66)];out[(an,b,'row31_sums')]=r31
    return out
f=feats(CT0)
for k,v in f.items(): print(k,v)
# null comparison for target hits and equal-diff runs
rng=random.Random(1);tg=0;run=[];chi=[]
for _ in range(2000):
    l=list(CT0);rng.shuffle(l);g=feats(''.join(l))
    tg+=sum(len(g[(a,b,'targets')]) for a in ('ST','KA') for b in (0,1));run.append(max(g[(a,b,'maxrun_equal_diff')] for a in ('ST','KA') for b in (0,1)));chi.append(max(g[(a,b,'mod10_chi')] for a in ('ST','KA') for b in (0,1)))
print('null: mean target-substring hits per text',tg/2000,'| K4:',sum(len(f[(a,b,'targets')]) for a in ('ST','KA') for b in (0,1)))
print('null: equal-diff run distribution max (median,95%)',np.median(run),np.percentile(run,95),'| K4:',max(f[(a,b,'maxrun_equal_diff')] for a in ('ST','KA') for b in (0,1)))
print('null: mod10 chi max (median,95%)',np.median(chi),np.percentile(chi,95),'| K4:',max(f[(a,b,'mod10_chi')] for a in ('ST','KA') for b in (0,1)))
