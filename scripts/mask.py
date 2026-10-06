import random
from collections import defaultdict
exec(open('board.py').read().split('# K4 cell')[0])
from stack import CT0,ANCH
KA=KRY;ST=STD;AN=sorted(ANCH)
# Kryptos tableau (right panel): header, 26 rows label+KA shifted (31 letters), footer
TAB=[' ABCDEFGHIJKLMNOPQRSTUVWXYZABCD']+[ST[r]+''.join(KA[(r+j)%26] for j in range(31)) for r in range(26)]+[' ABCDEFGHIJKLMNOPQRSTUVWXYZABCD']
PANEL=[l for l in LINES]
GRIDS={'tableau':TAB,'tableau_bottomup':TAB[::-1],'leftpanel':PANEL,'leftpanel_bottomup':PANEL[::-1]}
def keyneed(ct,A,rule):
    return {i:A[((A.index(ct[i])-A.index(ANCH[i])) if rule==0 else (A.index(ct[i])+A.index(ANCH[i])) if rule==1 else (A.index(ANCH[i])-A.index(ct[i])))%26] for i in AN}
def scan(ct,minred=6):
    hits=[]
    needs={(an,r):keyneed(ct,A,r) for an,A in (('KA',KA),('ST',ST)) for r in range(3)}
    for gname,G in GRIDS.items():
        R=len(G)
        for n in range(4,33):
            for o in range(n):
                rows_needed=(o+96)//n+1
                for t0 in range(0,R-rows_needed+1):
                    for serp in (0,1):
                        # slot of each anchor
                        slot={}
                        for i in AN:
                            t=(i+o)//n;s=(i+o)%n
                            if serp and t%2: s=n-1-s
                            slot[i]=(t,s)
                        groups=defaultdict(list)
                        for i in AN: groups[slot[i][1]].append(i)
                        red=len(AN)-len(groups)
                        if red<minred: continue
                        for key,need in needs.items():
                            ok=True;used=set()
                            for s,idx in groups.items():
                                W=min(len(G[t0+slot[i][0]]) for i in idx)
                                cols=[c for c in range(W) if all(G[t0+slot[i][0]][c]==need[i] for i in idx)]
                                if not cols: ok=False;break
                            if ok: hits.append((red,gname,n,o,t0,serp)+key)
    return hits
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    rng=random.Random(5);n=7;cols=sorted(rng.sample(range(1,32),n));t0=1
    key=[]
    for i in range(97):
        t=i//n;s=i%n
        if t%2: s=n-1-s
        key.append(TAB[t0+t][cols[s]])
    ctrl=''.join(KA[(KA.index(p)+KA.index(k))%26] for p,k in zip(pt,key))
    h=scan(ctrl);print('control (tableau, 7 holes/row, serpentine):',[x for x in h if x[1:6]==('tableau',7,0,1,1)][:2],'total hits',len(h),flush=True)
    nl=[]
    for sd in range(5):
        l=list(CT0);random.Random(sd).shuffle(l);nl.append(len(scan(''.join(l))))
    print('null hits per shuffle:',nl,flush=True)
    hk=scan(CT0);print('K4 hits',len(hk),sorted(hk,reverse=True)[:10])
