# Weltzeituhr mechanism as key source.
#  (A) the city text read face by face as the drum turns (running key), 2017 transcription + 1990-name variants
#  (B) hour-ring alignment keys / city-count keys are periodic-24 special cases (already excluded) -> checked directly anyway
import numpy as np,itertools,random,time,unicodedata
from stack import CT0,ANCH
from gen import colperms
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';AN=np.array(sorted(ANCH))
F={ -8:("VANCOUVER DAWSON SANFRANCISCO LOSANGELES",""), -7:("EDMONTON DENVER",""),
 -6:("NEWORLEANS MEXIKOSTADT","GUATEMALASTADT MANAGUA GALAPAGOSI"),
 -5:("MONTREAL WASHINGTON NEWYORK HAVANNA","PANAMA SANTAFEDEBOGOTA QUITO LIMA"),
 -4:("HALIFAX","CARACAS LAPAZ ASUNCION SANTIAGODECHILE"),
 -3:("WESTGRÖNLAND","BRASILIA RIODEJANEIRO SÃOPAULO MONTEVIDEO BUENOSAIRES"),
 -2:("OSTGRÖNLAND",""), -1:("AZOREN","KAPVERDE"),
 0:("REYKJAVIK DUBLIN LONDON LISSABON MADEIRA BISSAU","CASABLANCA CONAKRY DAKAR BAMAKO ACCRA"),
 1:("AMSTERDAM BERLIN BRÜSSEL BUDAPEST MADRID PARIS PRAG STOCKHOLM WARSCHAU","OSLO KOPENHAGEN WIEN BERN PRESSBURG BELGRAD ROM TUNIS KINSHASA"),
 2:("HELSINKI RIGA TALLINN WILNA MINSK KIEW BUKAREST SOFIA NIKOSIA","ANKARA ISTANBUL ATHEN TELAVIV JERUSALEM BEIRUT DAMASKUS KAIRO KAPSTADT"),
 3:("MURMANSK STPETERSBURG MOSKAU","TEHERAN BAGDAD ADEN SANAA ADDISABEBA MOGADISCHU DARESSALAM ANTANANARIVO KUWAIT"),
 4:("NISCHNIJNOWGOROD WOLGOGRAD BAKU TIFLIS ERIWAN","KABUL MAURITIUS"),
 5:("JEKATERINBURG ASCHGABAT BISCHKEK DUSCHANBE","NEWDELHI KARACHI COLOMBO"),
 6:("OMSK ALMATY TASCHKENT NOWOSIBIRSK","RANGUN DHAKA"), 7:("KRASNOJARSK","HANOI BANGKOK PHNOMPENH JAKARTA"),
 8:("IRKUTSK ULANBATOR","PEKING SHANGHAI MANILA PERTH HONGKONG KUALALUMPUR SINGAPUR"), 9:("JAKUTSK","PJÖNGJANG TOKYO SEOUL"),
 10:("CHABAROWSK WLADIWOSTOK","SYDNEY CANBERRA MELBOURNE"), 11:("MAGADAN SACHALIN",""), 12:("KAMTSCHATKA","DATUMSGRENZE WELLINGTON")}
GUESS={-9:("ANCHORAGE",""),-10:("HONOLULU",""),-11:("",""),-12:("","")}
OLD={'STPETERSBURG':'LENINGRAD','ALMATY':'ALMAATA','PRESSBURG':'BRATISLAVA','TELAVIV':'','JERUSALEM':'','NISCHNIJNOWGOROD':'GORKI','ULANBATOR':'ULANBATOR'}
def norm(s,mode):
    if mode=='oe': s=s.replace('Ö','OE').replace('Ü','UE').replace('Ä','AE')
    s=unicodedata.normalize('NFD',s);return ''.join(c for c in s if c.isalpha() and c.isascii()).upper()
def streams():
    out={}
    for era in ('2017','1990'):
        for g in (0,1):
            faces=dict(F);faces.update(GUESS if g else {})
            for mode in ('strip','oe'):
                for part in ('TB','BT','T','B'):
                    seq=[]
                    for z in range(-12,13):
                        if z not in faces: continue
                        t,b=faces[z]
                        if era=='1990':
                            t=' '.join(OLD.get(w,w) for w in t.split());b=' '.join(OLD.get(w,w) for w in b.split())
                        txt={'TB':t+' '+b,'BT':b+' '+t,'T':t,'B':b}[part]
                        seq.append(norm(txt,mode))
                    for d in ('east','west'):
                        s=''.join(seq if d=='east' else seq[::-1])
                        for wordrev in (0,1):
                            ss=s[::-1] if wordrev else s
                            if len(ss)>=30: out[(era,g,mode,part,d,wordrev)]=ss
    return out
def keys(ct):
    o={}
    for an,A in (('ST',ST),('KA',KA)):
        C=np.array([A.index(ct[i]) for i in AN]);P=np.array([A.index(ANCH[i]) for i in AN])
        o[(an,'vig')]=(C-P)%26;o[(an,'beau')]=(C+P)%26;o[(an,'var')]=(P-C)%26
    return o
def test(ct,S,perms=None,best_only=True):
    K=keys(ct);res=[]
    for lab,s in S.items():
        L=len(s)
        for (an,rule),k in K.items():
            A=ST if an=='ST' else KA
            v=np.array([A.index(c) for c in s]);vv=np.concatenate([v]*(int(200/L)+3))
            if perms is None:
                W=vv[np.arange(L)[:,None]+AN[None,:]];m=(W==k).sum(axis=1);o=int(np.argmax(m))
                res.append((int(m[o]),lab,an,rule,o))
            else:
                for c0 in range(0,len(perms),1500):
                    pos=perms[c0:c0+1500][:,AN].astype(np.int64)            # n x 24 ct positions
                    W=vv[np.arange(L)[None,:,None]+pos[:,None,:]]            # n x L x 24
                    m=(W==k[None,None,:]).sum(axis=2);x=np.unravel_index(np.argmax(m),m.shape)
                    res.append((int(m[x]),lab,an,rule,('perm',c0+int(x[0])),int(x[1])))
    res.sort(key=lambda r:-r[0]);return res
if __name__=='__main__':
    t=time.time();S=streams();print('streams',len(S),'lengths',sorted({len(v) for v in S.values()})[:3],'...',flush=True)
    SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    # control: running key = east TB 2017 strip stream from offset 37
    s=S[('2017',0,'strip','TB','east',0)];pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    ctl=''.join(KA[(KA.index(c)+KA.index(s[(i+37)%len(s)]))%26] for i,c in enumerate(pt))
    print('control direct',test(ctl,{('2017',0,'strip','TB','east',0):s})[:1],flush=True)
    for nm,ct in (('K4',CT0),('null',SH)):
        r=test(ct,S);print('(A) direct',nm,'best',r[:3],'| tests',len(r),flush=True)
    P=np.concatenate([colperms(97,W) for W in range(2,8)])
    S2={k:v for k,v in S.items() if k[2]=='strip' and k[5]==0 and k[3] in ('TB','T')}
    for nm,ct in (('K4',CT0),('null',SH)):
        r=test(ct,S2,P);print('(A) columnar 2-7',nm,'best',r[:3],round(time.time()-t),flush=True)
