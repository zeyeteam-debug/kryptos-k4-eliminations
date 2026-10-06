import random
exec(open('board.py').read().split('# K4 cell')[0])
from stack import CT0,ANCH
from sculpt import TAB
KA=KRY;ST=STD;AN=sorted(ANCH)
def need(ct,A,r):
    return {i:((A.index(ct[i])-A.index(ANCH[i])) if r==0 else (A.index(ct[i])+A.index(ANCH[i])) if r==1 else (A.index(ANCH[i])-A.index(ct[i])))%26 for i in AN}
def keystreams():
    # 1) YAR-shaped holes (relative cols 0,1,3) tiled every q columns, offset c0, over tableau / left panel, from any start row, rows top-down or bottom-up, straight or serpentine
    for gname,G0 in (('tableau',TAB),('panel',LINES)):
        for rev in (0,1):
            G=G0[::-1] if rev else G0
            for q in range(4,33):
                for c0 in range(q):
                    for t0 in range(len(G)):
                        for serp in (0,1):
                            ks=[];t=t0
                            while len(ks)<97 and t<len(G):
                                row=G[t];cols=[c for c in range(len(row)) if (c-c0)%q in (0,1,3) and row[c].isalpha()]
                                if serp and (t-t0)%2: cols=cols[::-1]
                                ks+=[row[c] for c in cols];t+=1
                            if len(ks)>=97: yield ('YARmask',gname,rev,q,c0,t0,serp),''.join(ks[:97])
    # 2) running key starting exactly at the raised Y (or A, R) and reading on (panel), and tableau from row N (same line as YAR)
    panel=''.join(LINES).replace('?','')
    y=panel.index('ENDYAHR')+3
    for s,nm in ((y,'fromY'),(y+1,'fromA'),(y+3,'fromR')):
        yield ('panel',nm),(panel[s:]+panel)[:97]
        yield ('panel_rev',nm),(panel[:s+1][::-1]+panel[::-1])[:97]
    tabtxt=''.join(TAB[14:]+TAB[:14]).replace(' ','')
    for c in range(32):
        yield ('tableau_fromN',c),(tabtxt[c:]+tabtxt)[:97]
def scan(ct):
    hits=[];N={(an,r):need(ct,A,r) for an,A in (('KA',KA),('ST',ST)) for r in range(3)};best=0
    for lab,ks in keystreams():
        for (an,r),nd in N.items():
            A=KA if an=='KA' else ST
            m=sum(A.index(ks[i])==nd[i] for i in AN);best=max(best,m)
            if m==24: hits.append((lab,an,r))
    return hits,best
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    for lab,ks in keystreams():
        if lab==('YARmask','tableau',0,9,2,1,1): break
    ctrl=''.join(KA[(KA.index(p)+KA.index(k))%26] for p,k in zip(pt,ks))
    h,b=scan(ctrl);print('control found:',[x for x in h if x[0]==lab][:1],'total',len(h))
    h,b=scan(CT0);print('K4 full matches:',h,'best partial match:',b,'/24')
    l=list(CT0);random.Random(3).shuffle(l);print('null best partial:',scan(''.join(l))[1],'/24')
