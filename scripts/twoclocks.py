# Coupling the two clocks: Mengenlehreuhr lamp state (24 lamps) selects Weltzeituhr faces (24 faces);
# the city text of the selected faces (in drum order) is the running key. Also the complementary (unlit) faces.
import numpy as np,random,time
from stack import CT0,ANCH
from clocks import mengen
from wcity import F,GUESS,OLD,norm,keys,ST,KA
AN=np.array(sorted(ANCH))
def facetext(z,era,part):
    if z not in F: return ''
    t,b=F[z]
    if era=='1990': t=' '.join(OLD.get(w,w) for w in t.split());b=' '.join(OLD.get(w,w) for w in b.split())
    return norm({'TB':t+' '+b,'T':t,'B':b}[part],'strip')
ZONES=list(range(-11,13))                    # 24 faces: -11..+12 (Datumsgrenze shares -12/+12)
def mappings():
    m={}
    for start in range(24):                  # lamp 0 sits on face ZONES[start]
        for d in (1,-1): m[(start,d)]=[ZONES[(start+d*j)%24] for j in range(24)]
    return m
MAPS=mappings();BER=ZONES.index(1)
def run(ct,starts=None,parts=('TB','T','B'),eras=('2017','1990')):
    K=keys(ct);best=[];top=[0];n=0;t=time.time()
    LUT={an:np.array([(ST if an=='ST' else KA).index(chr(65+i)) for i in range(26)]) for an in ('ST','KA')}
    cache={}
    for h in range(24):
        for mi in range(60):
            for s in (0,1):
                bits,_=mengen(h,mi,s)
                for (st,d),faces in MAPS.items():
                    if starts is not None and st not in starts: continue
                    for inv in (0,1):
                        sel=[z for b,z in zip(bits,faces) if b^inv]
                        for era in eras:
                            for part in parts:
                                key_=(tuple(sel),era,part)
                                if key_ in cache: continue
                                cache[key_]=1
                                txt=''.join(facetext(z,era,part) for z in sel)
                                if len(txt)<20: continue
                                L=len(txt);rep=np.concatenate([np.array([ord(c)-65 for c in txt])]*(200//L+2))
                                for (an,rule),k in K.items():
                                    v=LUT[an][rep];W=v[np.arange(L)[:,None]+AN[None,:]];m=(W==k).sum(axis=1);o=int(np.argmax(m));n+=L
                                    sc=int(m[o])
                                    if sc>=top[0]:
                                        if sc>top[0]: top[0]=sc;best.clear()
                                        if len(best)<5: best.append((sc,(h,mi,s),(st,d),inv,era,part,an,rule,o))
    return best,n
if __name__=='__main__':
    SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    # control: 14:37:s1, lamp0 on Berlin face going east, lit faces, 2017 TB, KA vig, offset 11
    bits,_=mengen(14,37,1);faces=MAPS[(BER,1)];sel=[z for b,z in zip(bits,faces) if b]
    txt=''.join(facetext(z,'2017','TB') for z in sel)
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    ctl=''.join(KA[(KA.index(c)+KA.index(txt[(i+11)%len(txt)]))%26] for i,c in enumerate(pt))
    b,_=run(ctl,starts=[BER],parts=('TB',),eras=('2017',));print('control',b[0],flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        b,n=run(ct);print(nm,'window tests',n,'best',b[:4],flush=True)
