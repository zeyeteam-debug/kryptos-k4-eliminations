import numpy as np,json,itertools,sys,random
from stack import CT0,ANCH
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
WORDS=sorted(set('''KRYPTOS PALIMPSEST ABSCISSA ORDINATE BERLIN CLOCK BERLINCLOCK EASTNORTHEAST NORTHEAST EAST SHADOW SHADOWS LIGHT
LODESTONE MAGNETIC MAGNETISM COMPASS WEBSTER WILLIAM SANBORN SCHEIDT LANGLEY CIA LAYER LAYERTWO LAYERS IQLUSION ILLUSION UNDERGRUUND UNDERGROUND
DESPARATLY DESPERATELY CARTER TUTANKHAMUN TOMB CANDLE CHAMBER EGYPT NUANCE SUBTLE ABSENCE INVISIBLE LOCATION BURIED MESSAGE
WELTZEITUHR URANIA ALEXANDERPLATZ MOSKAU MOSCOW WALL TIME HOUR POSITION COORDINATE SECRET HIDDEN KEY CIPHER MATRIX MASK MASKING
DIGETAL INTERPRETATIU VIRTUALLY LUCID MEMORY FORCES YAR RAY HILL DELTA POINT ORIGIN QUADRANT ENIGMA PUZZLE SCULPTURE COPPER
PETRIFIED WOOD QUARTZ GRANITE WATER POOL JIM JAMES ED EDWARD'''.split()))
CONTROL_WORDS=['UNDERLAY','OCHRE','VERDIGRIS']
def order_of(key): return sorted(range(len(key)),key=lambda i:(key[i],i))
def colidx(n,key):
    w=len(key);idx=[]
    for c in order_of(key): idx+=list(range(c,n,w))
    return np.array(idx)          # out[k]=in[idx[k]]
def search(ct,anch,words,AL=KA,label=''):
    n=len(ct);AN=np.array(sorted(anch));P=np.array([AL.index(anch[i]) for i in AN]);C=np.array([AL.index(c) for c in ct])
    W=[(w,np.array([AL.index(c) for c in w])) for w in words]
    hits=[]
    # transposition options: none, T(w), inverse T(w)
    Topts=[('none',np.arange(n))]
    for w,_ in W:
        if len(w)>1:
            idx=colidx(n,w);inv=np.argsort(idx)
            Topts.append(('T:'+w,inv))           # pos of pt i in ct = inv[i] (T applied to pt first)
            Topts.append(('Tinv:'+w,idx))
    # keys of words at given positions
    for tname,pos in Topts:
        q=pos[AN]                                 # ct positions of anchors
        cq=C[q]
        for s1 in (1,-1):
            # needed combined key (in ct-index frame, since Q applied after T): c = p + s1*k1 + s2*k2
            need=(cq-P)%26
            K1=np.stack([wv[q%len(wv)] for _,wv in W])          # nW x 24
            for s2 in (1,-1):
                # need - s1*K1  must equal s2*K2
                lhs=(need[None,:]-s1*K1)%26                     # nW x 24
                rhs=(s2*K1)%26                                  # nW x 24 (K2 from same list)
                # match rows: use hashing
                d={}
                for j in range(len(W)): d.setdefault(rhs[j].tobytes(),[]).append(j)
                for i in range(len(W)):
                    for j in d.get(lhs[i].tobytes(),[]):
                        hits.append((tname,W[i][0],s1,W[j][0],s2))
    return hits
def search_TTQ(ct,anch,words,AL=KA):
    """PK6 style: T(w1) -> T(w2) -> Q(w3)"""
    n=len(ct);AN=np.array(sorted(anch));P=np.array([AL.index(anch[i]) for i in AN]);C=np.array([AL.index(c) for c in ct])
    W=[(w,np.array([AL.index(c) for c in w])) for w in words if len(w)>1]
    inv={w:np.argsort(colidx(n,w)) for w,_ in W}
    hits=[]
    for w1,_ in W:
        for w2,_ in W:
            q=inv[w2][inv[w1][AN]]
            need=(C[q]-P)%26
            for w3,v3 in W:
                for s in (1,-1):
                    if np.array_equal((s*v3[q%len(v3)])%26,need): hits.append((w1,w2,w3,s))
    return hits
if __name__=='__main__':
    d=json.load(open('pk_controls.json'))
    pk4=open('pk4.txt').read().replace('\n','').strip();p4=d['PK4']['plaintext']
    anch4={i:p4[i] for i in list(range(21,34))+list(range(63,74))}
    words=WORDS+CONTROL_WORDS
    print('control PK4 (24 anchors):',[h for h in search(pk4,anch4,words) if 'UNDERLAY' in h[0]][:3],flush=True)
    pk6=open('pk6.txt').read().replace('\n','').strip();p6=d['PK6']['plaintext']
    anch6={i:p6[i] for i in list(range(21,34))+list(range(63,74))}
    print('control PK6 (24 anchors):',[h for h in search_TTQ(pk6,anch6,WORDS+['HANDIWORK','SMITHWORK','PORTAL']) if h[2]=='PORTAL'][:2],flush=True)
    print('words',len(WORDS))
    for AL,an in ((KA,'KA'),(ST,'ST')):
        h=search(CT0,ANCH,WORDS,AL);print(an,'K4 T+Q+Q hits:',len(h),h[:10],flush=True)
        h=search_TTQ(CT0,ANCH,WORDS,AL);print(an,'K4 T+T+Q hits:',len(h),h[:10],flush=True)
