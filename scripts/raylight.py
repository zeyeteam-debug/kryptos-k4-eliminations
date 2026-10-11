# "Look through the YAR hole toward ENE": read K4 along parallel light rays of a compass bearing.
# K4 laid in its carved 4-line layout (from the ? on line 25) or in rows of width w; rays of slope dy/dx;
# order = ray by ray (sorted by intercept), along the ray, both senses; every compass point of 16 + ENE approximations.
exec(open('handmask.py').read().split('# letter cells')[0].replace("print(","(lambda *a,**k:None)("))
import random
DIRS={'E':(1,0),'ENE_1:2':(2,1),'ENE_2:5':(5,2),'ENE_5:12':(12,5),'ENE_1:3':(3,1),'NE':(1,1),'NNE':(1,2),'NNE_2:5':(2,5),'N':(0,1),
      'NNW':(-1,2),'NW':(-1,1),'WNW':(-2,1),'WNW_5:12':(-12,5),'ESE':(2,-1),'SE':(1,-1),'SSE':(1,-2)}
def layouts():
    yield 'carved',[(L-24,c) for L,c in cells]
    for w in range(4,49): yield f'w{w}',[(i//w,i%w) for i in range(N)]
    for w in range(4,49):
        for o in range(1,w): yield f'w{w}o{o}',[((i+o)//w,(i+o)%w) for i in range(N)]
P=[];Lb=[]
for ln,pos in layouts():
    for dn,(dx,dy) in DIRS.items():
        # 'up' = north = decreasing row; ray parameter t along (dx,-dy); intercept = cross product
        key_line=[pos[i][1]*dy+pos[i][0]*dx for i in range(N)]   # constant along the ray
        key_t=[pos[i][1]*dx-pos[i][0]*dy for i in range(N)]
        for s1 in (1,-1):
            for s2 in (1,-1):
                o=sorted(range(N),key=lambda i:(s1*key_line[i],s2*key_t[i],i))
                o=np.array(o);P+=[o,np.argsort(o)];Lb+=[(ln,dn,s1,s2,'a'),(ln,dn,s1,s2,'b')]
P=np.array(P,dtype=np.int16);u,ix=np.unique(P,axis=0,return_index=True);ix=sorted(ix);P=P[ix];Lb=[Lb[i] for i in ix]
print('ray perms',len(P),flush=True)
SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
for a,ch in ANCH.items(): pt[a]=ch
i=[x for x,l in enumerate(Lb) if l[0]=='carved' and l[1].startswith('ENE')][0];m=P[i];ct=['?']*N
for a in range(N): ct[m[a]]=KA[(KA.index(pt[a])+KA.index('LIGHTRAY'[a%8]))%26]
print('control',[h for h in run(''.join(ct),ANCH,P[i:i+1],[Lb[i]],Pmax=12) if h[0]=='key'][:1],flush=True)
for nm,x in (('K4',CT0),('null',SH)):
    H=[];qq=0
    for c0 in range(0,len(P),20000):
        H+=run(x,ANCH,P[c0:c0+20000],Lb[c0:c0+20000],Pmax=24);qq+=qii(x,P[c0:c0+20000],24)
    print(nm,'periodic',len(H),[h for h in H if h[0]=='key'][:3],'any',sum(1 for h in H if h[0]=='any'),'QI/II',qq,flush=True)
