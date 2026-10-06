import numpy as np,random
from stack import CT0,ANCH
from gen import run
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';AN=sorted(ANCH);N=97
MORSE_TXT='VIRTUALLYINVISIBLEDIGETALINTERPRETATIUSHADOWFORCESLUCIDMEMORYTISYOURPOSITIONSOSRQ'
MC={'A':'.-','B':'-...','C':'-.-.','D':'-..','E':'.','F':'..-.','G':'--.','H':'....','I':'..','J':'.---','K':'-.-','L':'.-..','M':'--','N':'-.','O':'---','P':'.--.','Q':'--.-','R':'.-.','S':'...','T':'-','U':'..-','V':'...-','W':'.--','X':'-..-','Y':'-.--','Z':'--..'}
codes=[MC[c] for c in MORSE_TXT]
seqs={'morse_len':[len(c) for c in codes],'morse_dots':[c.count('.') for c in codes],'morse_dashes':[c.count('-') for c in codes],
 'morse_binval':[int(c.replace('.','0').replace('-','1'),2)+ (1<<len(c)) for c in codes]}
bits=''.join(c.replace('.','0').replace('-','1') for c in codes)
seqs['morse_bits5']=[int(bits[i:i+5],2) for i in range(0,len(bits)-4,5)]
seqs['morse_bits_sliding5']=[int(bits[i:i+5],2) for i in range(len(bits)-4)]
# clock: offsets of cities in face order (Europe/Asia faces we transcribed), one number per city, and per letter
faces={-1:'AZOREN KAPVERDE',0:'REYKJAVIK DUBLIN LONDON LISSABON MADEIRA BISSAU CASABLANCA CONAKRY DAKAR BAMAKO ACCRA',
1:'AMSTERDAM BERLIN BRUSSEL BUDAPEST MADRID PARIS PRAG STOCKHOLM WARSCHAU OSLO KOPENHAGEN WIEN BERN PRESSBURG BELGRAD ROM TUNIS KINSHASA',
3:'MOSKAU',5:'JEKATERINBURG BISCHKEK DUSCHANBE NEWDELHI KARACHI COLOMBO',6:'OMSK ALMATY TASCHKENT NOWOSIBIRSK RANGUN DHAKA',7:'KRASNOJARSK HANOI BANGKOK PHNOMPENH JAKARTA',
8:'IRKUTSK ULANBATOR PEKING SHANGHAI MANILA PERTH HONGKONG KUALALUMPUR SINGAPUR',9:'JAKUTSK PJONGJANG TOKYO SEOUL',10:'CHABAROWSK WLADIWOSTOK SYDNEY CANBERRA MELBOURNE',11:'MAGADAN SACHALIN',12:'KAMTSCHATKA DATUMSGRENZE WELLINGTON'}
seqs['clock_city_offsets']=[z for z,f in faces.items() for _ in f.split()]
seqs['clock_letter_offsets']=[z for z,f in faces.items() for w in f.split() for _ in w]
seqs['clock_namelens']=[len(w) for z,f in faces.items() for w in f.split()]
seqs['Wseg_lengths']=[20,15,11,9,15,22];seqs['Wpos']=[21,37,49,59,75]
def scan(ct):
    hits=[];best=0;n=0
    for nm,s0 in seqs.items():
        for rev in (0,1):
            s=s0[::-1] if rev else s0
            ext=s*(97//len(s)+3)
            for mul in (1,-1,2,3):
                for an,A in (('ST',ST),('KA',KA)):
                    c=np.array([A.index(ct[i]) for i in AN]);p=np.array([A.index(ANCH[i]) for i in AN])
                    for rule in (0,1):
                        for off in range(len(s)):
                            k=np.array([(mul*ext[off+i])%26 for i in AN])
                            exp=(p+k)%26 if rule==0 else (k-p)%26
                            m=int((exp==c).sum());best=max(best,m);n+=1
                            if m==24: hits.append((nm,rev,mul,an,rule,off))
    return hits,best,n
# irregular rows from W segment lengths -> read by columns (transposition) + periodic keys
def rowsT(lengths,colorder=None,updown=0):
    rows=[];s=0
    for L in lengths: rows.append(list(range(s,min(N,s+L))));s+=L
    while s<N: rows.append(list(range(s,min(N,s+lengths[-1]))));s+=lengths[-1]
    W=max(len(r) for r in rows);out=[]
    for c in (colorder or range(W)):
        col=[r[c] for r in rows if c<len(r)]
        out+=col[::-1] if updown else col
    m=np.empty(N,dtype=np.int16);m[np.array(out)]=np.arange(N);return m
if __name__=='__main__':
    src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
    pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
    s=seqs['morse_len'];ext=s*3
    ctrl=''.join(KA[(KA.index(x)+ext[10+i])%26] for i,x in enumerate(pt))
    print('control morse_len:',scan(ctrl)[0][:1])
    h,b,n=scan(CT0);print('K4 numeric keys: hits',h,'best',b,'/24 of',n)
    l=list(CT0);random.Random(1).shuffle(l);print('null best',scan(''.join(l))[1])
    Ms=[];labs=[]
    for lens in ([20,15,11,9,15,22],[20,1,15,1,11,1,9,1,15,1,22],[21,16,12,10,16,22]):
        for ud in (0,1):
            m=rowsT(lens,updown=ud);Ms+=[m,np.argsort(m).astype(np.int16)];labs+=[(tuple(lens),ud,'enc'),(tuple(lens),ud,'dec')]
    print('W-row transpositions + any periodic key:',run(CT0,ANCH,np.array(Ms),labs,Pmax=12,RTHR=5,ANYSTRONG=3))
