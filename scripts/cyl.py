import numpy as np,random,sys,itertools
from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
N=97;AN=np.array(sorted(ANCH))
def colseq(W,s,k):
    seen=[];start=s
    while len(seen)<W:
        c=start
        while c not in seen: seen.append(c);c=(c+k)%W
        start=(start+1)%W
        while start in seen and len(seen)<W: start=(start+1)%W
    return seen
def perm(W,s,k,h,cd):
    cells={}
    for i in range(N):
        r=i//W;c=(i%W+r*h)%W;cells.setdefault(c,[]).append((r,i))
    out=[]
    for n,c in enumerate(colseq(W,s,k)):
        col=sorted(cells.get(c,[]))
        if cd==1 or (cd==2 and n%2): col=col[::-1]
        out+=[i for _,i in col]
    m=np.empty(N,dtype=np.int16);m[np.array(out)]=np.arange(N);return m   # m[pt]=ct pos
def family():
    labs=[];Ms=[]
    for W in (24,12):
        for s in range(W):
            for k in range(1,W):
                for h in range(W):
                    for cd in (0,1,2):
                        m=perm(W,s,k,h,cd);Ms.append(m);labs.append((W,s,k,h,cd,'enc'))
                        Ms.append(np.argsort(m).astype(np.int16));labs.append((W,s,k,h,cd,'dec'))
    return np.array(Ms),labs
RTHR=6
def run(ct,Ms,labs):
    hits=[];pairs=np.array(list(itertools.combinations(range(len(AN)),2)))
    for an,A in (('ST',ST),('KA',KA)):
        ci=np.array([A.index(x) for x in ct]);P=np.array([A.index(ANCH[i]) for i in AN])
        C=ci[Ms[:,AN]]
        for rule in (0,1):
            K=(C-P)%26 if rule==0 else (C+P)%26
            same=K[:,pairs[:,0]]==K[:,pairs[:,1]]
            for Per in range(1,13):
                # key indexed by PT position
                res=AN%Per;R=len(AN)-len(set(res))
                if R>=RTHR:
                    need=res[pairs[:,0]]==res[pairs[:,1]]
                    ok=np.all(same[:,need],axis=1)
                    for x in np.nonzero(ok)[0]: hits.append(('ptkey',Per,an,rule,labs[x]))
                # key indexed by CT position
                rs=Ms[:,AN]%Per
                need=rs[:,pairs[:,0]]==rs[:,pairs[:,1]]
                ok=np.all(same|~need,axis=1)
                srt=np.sort(rs,axis=1);distinct=1+np.sum(srt[:,1:]!=srt[:,:-1],axis=1)
                ok&=(len(AN)-distinct)>=RTHR
                for x in np.nonzero(ok)[0]: hits.append(('ctkey',Per,an,rule,labs[x]))
    return hits
if __name__=='__main__':
    Ms,labs=family();print('permutations',len(Ms),flush=True)
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    m=perm(24,1,5,3,2);key='KRYPTOS';ct=['?']*N
    for i in range(N): ct[m[i]]=KA[(KA.index(pt[i])+KA.index(key[i%7]))%26]
    h=run(''.join(ct),Ms,labs);print('control found:',[x for x in h if x[4]==(24,1,5,3,2,'enc')][:2],'total hits',len(h))
    random.seed(3);s=list(CT0);random.shuffle(s);print('null hits:',len(run(''.join(s),Ms,labs)))
    h=run(CT0,Ms,labs);print('K4 hits',len(h),h[:20])
