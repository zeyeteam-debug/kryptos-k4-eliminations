import numpy as np
exec(open('board.py').read().split('# K4 cell')[0])
allct=''.join(LINES)
K3C=allct[allct.index('ENDYAH'):allct.index('?OBKR')].replace('?','')
def rot(N,w,cw=True):
    # returns idx: out[n]=in[idx[n]]  ; write rows width w (ragged), rotate, read rows
    rows=(N+w-1)//w;idx=[]
    cols=range(w) if cw else range(w-1,-1,-1)
    for c in cols:
        rr=range(rows-1,-1,-1) if cw else range(rows)
        for r in rr:
            i=r*w+c
            if i<N: idx.append(i)
    return np.array(idx)
def compose(N,w1,d1,w2,d2):
    a=rot(N,w1,d1);b=rot(N,w2,d2);return a[b]     # out[n]=in[a[b[n]]]
if __name__=='__main__':
    print(len(K3C),len(K3P))
    N=len(K3C);P=np.array(list(K3P[:N]));found=[]
    for w1 in range(2,60):
        for w2 in range(2,60):
            for d1 in (True,False):
                for d2 in (True,False):
                    idx=compose(N,w1,d1,w2,d2)
                    s=''.join(P[idx])
                    m=sum(a==b for a,b in zip(s,K3C))
                    if m>N*0.9: found.append((m,w1,d1,w2,d2))
    print(sorted(found,reverse=True)[:5])
