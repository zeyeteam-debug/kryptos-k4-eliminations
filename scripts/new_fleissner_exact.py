# Exact (exhaustive) search: 10x10 Fleissner turning grille x periodic substitution, K4 cribs.
# Conventions: page frame fixed; grille hole set H0 (one cell per 90-deg orbit); rotation t holes H_t=R^t(H0)
# (R clockwise (r,c)->(c,9-r); ccw uses R^-1). Plaintext 0..99 written: t=0..3, holes of H_t in WRITE order
# (row-major or col-major of page). Blanks = pt 97,98,99 (end) or pt 0,1,2 (start). Ciphertext = filled cells
# in READ order (row/col-major), blanks skipped. Key indexed by K4 plaintext index i (0..96).
# Models: {ST,KA}x{Vig c=p+k, Beau c=k-p} periods 1..12 (48), + arbitrary-alphabet-per-class periods in ARB.
import subprocess,random,sys,time,os
BIN=os.environ.get('FLBIN','./fl_exact7')
# model masks (hex): bits a*24+t*12+(p-1) for Vig/Beau, bit 48 = ARB-1
MASK_P8=hex(sum(1<<(a*24+t*12+p-1) for a in (0,1) for t in (0,1) for p in range(1,9))|(1<<48))[2:]
from stack import CT0,ANCH
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
ROW=list(range(100));COL=[r*10+c for c in range(10) for r in range(10)]
def R(x,d):
    r,c=divmod(x,10)
    return c*10+9-r if d==1 else (9-c)*10+r
def model_names(arb):
    n=[]
    for m in range(48):
        a,t,p=m//24,(m//12)%2,m%12+1;n.append('%s-%s-%d'%('ST' if a==0 else 'KA','V' if t==0 else 'B',p))
    return n+['ARB-%d'%p for p in arb]
def encrypt(H0,pt100,d,W,RDo):
    wp={x:i for i,x in enumerate(W)};cell={};H=list(H0);k=0
    for t in range(4):
        for x in sorted(H,key=lambda x:wp[x]): cell[x]=pt100[k];k+=1
        H=[R(x,d) for x in H]
    assert len(cell)==100
    return ''.join(cell[x] for x in RDo if cell[x]!='.')
def run(ct,d,W,RDo,bmode,arb=(1,),mask=None):
    inp='%s\n%d %d %d %s\n%s\n%s\n'%(ct,d,bmode,len(arb),' '.join(map(str,arb)),' '.join(map(str,W)),' '.join(map(str,RDo)))
    t=time.time()
    r=subprocess.run(['nice','-n','5',BIN]+([mask] if mask else []),input=inp,capture_output=True,text=True)
    sols=[];names=model_names(arb)
    for line in r.stdout.splitlines():
        if line.startswith('SOL'):
            al=int(line.split()[1][6:],16);H0=[int(v) for v in line.split('H0=')[1].strip(',').split(',')]
            sols.append((H0,[names[m] for m in range(len(names)) if al>>m&1],0,None))
        elif line.startswith('FAM'):
            al=int(line.split()[1][6:],16);fr=int(line.split()[2][5:]);lab=line.split('lab=')[1].strip()
            sols.append(([x for x in range(100) if lab[x]=='0'],[names[m] for m in range(len(names)) if al>>m&1],fr,lab))
        elif line.startswith('DONE'): done=line
    return sols,done,time.time()-t
def planted(sol,H0):
    # family contains planted grille iff all determined cells with label 0 match and no planted hole is labelled otherwise
    H0s=set(H0)
    if sol[3] is None: return set(sol[0])==H0s
    return all((sol[3][x]=='0')==(x in H0s) for x in range(100) if sol[3][x]!='.')
VARIANTS=[(d,wn,rn,b) for b in (0,1) for d in (1,-1) for wn,rn in (('row','row'),('row','col'))]
OD={'row':ROW,'col':COL}
if __name__=='__main__':
    what=sys.argv[1]
    if what=='ctrl':
        # planted control: random grille, random plaintext with cribs, KA Vigenere period 5
        seed=int(sys.argv[2]);vi=int(sys.argv[3]);d,wn,rn,b=VARIANTS[vi];mask=sys.argv[4] if len(sys.argv)>4 else None
        rng=random.Random(seed)
        orbs=[];seen=set()
        for x in range(100):
            if x in seen: continue
            o=[x]
            for _ in range(3): o.append(R(o[-1],1))
            seen.update(o);orbs.append(o)
        H0=sorted(rng.choice(o) for o in orbs)
        k4=[rng.choice(ST) for _ in range(97)]
        for i,c in ANCH.items(): k4[i]=c
        key=[rng.randrange(26) for _ in range(5)]
        ctk=[KA[(KA.index(c)+key[i%5])%26] for i,c in enumerate(k4)]
        pt100=ctk+['.']*3 if b==0 else ['.']*3+ctk
        ct=encrypt(H0,pt100,d,OD[wn],OD[rn])
        sols,done,dt=run(ct,d,OD[wn],OD[rn],b,mask=mask)
        print('CONTROL seed',seed,'variant dir=%s write=%s read=%s blanks=%s'%('cw' if d==1 else 'ccw',wn,rn,'end' if b==0 else 'start'),'mask',mask)
        print('  planted H0',H0,'key',key);print(' ',done.strip(),'time %.1fs'%dt)
        for s_ in sols: print('  survivor',planted(s_,H0) and 'PLANTED' or 'other',s_[1],'free orbits',s_[2])
        print('  planted found with KA-V-5 alive:',any(planted(s_,H0) and 'KA-V-5' in s_[1] for s_ in sols),flush=True)
    else:
        if what=='K4': ct=CT0
        else:
            l=list(CT0);random.Random(int(what[4:])).shuffle(l);ct=''.join(l)
        vs=VARIANTS if len(sys.argv)<3 else [VARIANTS[int(v)] for v in sys.argv[2].split(',')]
        mask={'all':None,'p8':MASK_P8}[sys.argv[3]] if len(sys.argv)>3 else None
        for d,wn,rn,b in vs:
            sols,done,dt=run(ct,d,OD[wn],OD[rn],b,mask=mask)
            print(what,'models=%s'%(sys.argv[3] if len(sys.argv)>3 else 'all'),'dir=%s write=%s read=%s blanks=%s'%('cw' if d==1 else 'ccw',wn,rn,'end' if b==0 else 'start'),done,'%.0fs'%dt,'survivors',len(sols),flush=True)
            for s in sols: print('    H0(partial)=',s[0],'free_orbits=',s[2],'models=',s[1],'labels=',s[3],flush=True)
