import numpy as np,random,string,itertools
from stack import CT0,ANCH
from gen import run
base=['TRAIN','TRAINS','RAIL','RAILS','RAILWAY','RAILROAD','TRACK','TRACKS','LINE','LINES','STATION','PLATFORM','TRAM','SBAHN','UBAHN','BAHNHOF',
      'TOWARDS','TOWARD','POINTING','POINTS','POINTSTO','LEADS','LEADING','DIRECTION','HEADING','FACING','FROM','TO','OF','AT','NEAR','PAST']
arts=['','THE','A','AN','DER','DIE','DAS']
W=sorted({a+b for a in arts for b in base})
def ok_any(an,P):
    return [x for x in run(CT0,an,np.arange(97,dtype=np.int16)[None,:],['g'],Pmax=24,RTHR=4,ANYSTRONG=2)]
cnt=0;hits=[];hk=[]
for w in W:
    for s in range(34,63-len(w)+1):
        an=dict(ANCH);an.update({s+j:c for j,c in enumerate(w)});cnt+=1
        H=ok_any(an,24)
        for x in H:
            if x[0]=='key': hk.append((w,s,x))
            elif x[2]==19: hits.append((w,s,x[3]))
# random control of same lengths
rng=random.Random(3);rc=0;rh=0
for w in W:
    r=''.join(rng.choice(string.ascii_uppercase) for _ in w)
    for s in range(34,63-len(r)+1):
        an=dict(ANCH);an.update({s+j:c for j,c in enumerate(r)});rc+=1
        rh+=any(x[0]=='any' and x[2]==19 for x in ok_any(an,24))
print('gap fillers',len(W),'placements',cnt,'ST/KA periodic (<=24) consistent:',len(hk),hk[:3])
print('period-19 any-alphabet consistent:',len(hits),'rate',round(len(hits)/cnt,3),'| random strings rate',round(rh/rc,3))
