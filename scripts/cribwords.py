import numpy as np,random,time,itertools,string
from stack import CT0,ANCH
from gen import run
from greek import ALS
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
# (a) the clue words themselves as keys in all PK-style pipelines
CW=['EASTNORTHEAST','BERLINCLOCK','BERLIN','CLOCK','EAST','NORTHEAST','ENE','EASTNORTHEASTBERLINCLOCK','BERLINCLOCKEASTNORTHEAST','NORTH','EASTBERLIN','WELTZEITUHR','ALEXANDERPLATZ']
# (b) Berlin / ENE-bearing words as extra cribs at every free position; test period-19 (and 1-24) with ST/KA and arbitrary alphabets
BW=['LENIN','LENINPLATZ','LENINALLEE','MOLLSTRASSE','ALEXANDERPLATZ','ALEXANDER','FERNSEHTURM','WELTZEITUHR','URANIA','KARLMARXALLEE',
    'BAHNHOF','SBAHN','UBAHN','STATION','TRAIN','TRAM','LINE','MAUER','WALL','THEWALL','BERLINWALL','CHECKPOINT','CHARLIE','CHECKPOINTCHARLIE',
    'BRANDENBURG','GATE','TOR','BRANDENBURGERTOR','UNTERDENLINDEN','REICHSTAG','TOWER','TVTOWER','WORLDCLOCK','HOUSE','HAUS','HAUSDESLEHRERS',
    'HAUSDESREISENS','DEGREES','SIXTYSEVEN','BEARING','COMPASS','SHADOW','LIGHT','TIME','HOUR','MINUTE','MIDNIGHT','NOVEMBER','NINTH',
    'NINETEENEIGHTYNINE','FALL','DIVIDED','EAST','WEST','EASTERN','WESTERN','GERMANY','DDR','GDR','SOVIET','KGB','STASI','AGENT','LAYER','TWO','THREE']
AN0=set(ANCH)
def free_slots(L):
    return [s for s in range(0,98-L) if not (set(range(s,s+L))&AN0)]
def guess_test(ct,words,Pset=(19,)):
    out={}
    for w in words:
        n=0;h=0;best=None
        for s in free_slots(len(w)):
            an=dict(ANCH);an.update({s+j:c for j,c in enumerate(w)})
            H=run(ct,an,np.arange(97,dtype=np.int16)[None,:],[(w,s)],Pmax=24,RTHR=6,ANYSTRONG=6)
            H=[x for x in H if (x[0]=='key' and x[4] in Pset) or (x[0]=='any' and x[2] in Pset)];n+=1
            if H: h+=1;best=H[0]
        out[w]=(n,h,best)
    return out
if __name__=='__main__':
    t=time.time()
    r=guess_test(CT0,BW);tot=sum(v[0] for v in r.values());hits={k:v for k,v in r.items() if v[1]}
    print('(b) K4 placements',tot,'consistent with period 19 (>=6 checks):',hits,round(time.time()-t),flush=True)
    rng=random.Random(5);RW=[''.join(rng.choice(string.ascii_uppercase) for _ in range(len(w))) for w in BW for _ in range(3)]
    rr=guess_test(CT0,RW);print('(b) random-word control: placements',sum(v[0] for v in rr.values()),'consistent',sum(v[1] for v in rr.values()),flush=True)
    # also all periods 1-24 for the Berlin words
    r2=guess_test(CT0,BW,Pset=tuple(range(1,25)));print('(b2) any period<=24:',{k:v for k,v in r2.items() if v[1]},flush=True)
    import dictsearch as D
    routes=[('none',np.arange(97))]
    for Wd in range(2,49):
        for corner in ('TL','TR','BL','BR'):
            for cw in (0,1):
                m=D.perm_from_read(D.spiral_order(97,Wd,corner,cw),97);routes+=[((Wd,corner,cw,'r'),m),((Wd,corner,cw,'w'),np.argsort(m))]
    for an,AL in ALS.items():
        TB=D.table(CW+BW,ANCH,AL);h=D.run(CT0,ANCH,routes,[2,3,4,5,6,7],AL,TB);print('(a) Q->spiral->col',an,len(h),h[:3],flush=True)
    from layers3 import search,search_TTQ
    from layers3b import search_QQT
    for an,AL in ALS.items():
        W=CW+BW;h=search(CT0,ANCH,W,AL);h2,h3=search_QQT(CT0,ANCH,W,AL);h4=search_TTQ(CT0,ANCH,W,AL)
        print('(a) layers',an,len(h),len(h2),len(h3),len(h4),(h+h2+h3+h4)[:3],round(time.time()-t),flush=True)
