import numpy as np,itertools
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
def run(ct,anch,Ms,labs,Pmax=12,RTHR=6,anyalpha=True,ANYSTRONG=4):
    """Ms[n][pt]=ct index. anchors: dict pt_index->letter. '?' in ct = unknown (anchor dropped)."""
    AN=np.array(sorted(i for i in anch if True));hits=[]
    if len(AN)<2: return hits
    PR=np.array(list(itertools.combinations(range(len(AN)),2)))
    CI=Ms[:,AN];cl=np.array([ord(x) for x in ct]);Cl=cl[CI];valid=Cl!=63   # '?'
    vpair=valid[:,PR[:,0]]&valid[:,PR[:,1]]
    for an,A in (('ST',ST),('KA',KA)):
        lut=np.full(128,0);lut[[ord(c) for c in A]]=np.arange(26)
        C=lut[Cl];P=np.array([A.index(anch[i]) for i in AN])
        for rule in (0,1):
            K=(C-P)%26 if rule==0 else (C+P)%26
            same=K[:,PR[:,0]]==K[:,PR[:,1]]
            for Per in range(1,Pmax+1):
                for kind,R in (('pt',np.broadcast_to(AN%Per,CI.shape)),('ct',CI%Per)):
                    need=(R[:,PR[:,0]]==R[:,PR[:,1]])&vpair
                    ok=np.all(same|~need,axis=1)
                    # redundancy = valid anchors - distinct residues among valid
                    if ok.any():
                        for x in np.nonzero(ok)[0]:
                            rv=R[x][valid[x]];red=len(rv)-len(set(rv.tolist()))
                            if red>=RTHR: hits.append(('key',an,rule,kind,Per,red,labs[x]))
    if anyalpha:
        Pl=np.array([ord(anch[i]) for i in AN]);Peq=Pl[PR[:,0]]==Pl[PR[:,1]]
        Ceq=Cl[:,PR[:,0]]==Cl[:,PR[:,1]];bad=(Peq^Ceq)&vpair
        for Per in range(1,Pmax+1):
            for kind,R in (('pt',np.broadcast_to(AN%Per,CI.shape)),('ct',CI%Per)):
                sm=(R[:,PR[:,0]]==R[:,PR[:,1]])
                ok=~np.any(bad&sm,axis=1)&((sm&Peq&vpair).sum(axis=1)>=ANYSTRONG)
                for x in np.nonzero(ok)[0]: hits.append(('any',kind,Per,int((sm[x]&Peq&vpair[x]).sum()),labs[x]))
    return hits
def colperms(N,W):
    out=[]
    L=N//W;r=N%W;lens=[L+(1 if c<r else 0) for c in range(W)]
    for o in itertools.permutations(range(W)):
        idx=[]
        for c in o: idx+=[row*W+c for row in range(lens[c])]
        m=np.empty(N,dtype=np.int16);m[np.array(idx)]=np.arange(N);out.append(m)
    M=np.array(out);return np.concatenate([M,np.argsort(M,axis=1).astype(np.int16)])
