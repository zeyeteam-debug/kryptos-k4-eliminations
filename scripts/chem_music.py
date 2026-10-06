import numpy as np,random
from stack import CT0,ANCH
from gen import colperms
from rk_trans import prep,search
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';AN=sorted(ANCH)
SYM='H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og'.split()
MASS=[1,4,7,9,11,12,14,16,19,20,23,24,27,28,31,32,35,40,39,40,45,48,51,52,55,56,59,59,64,65,70,73,75,79,80,84,85,88,89,91,93,96,98,101,103,106,108,112,115,119,122,128,127,131,133,137,139,140,141,144,145,150,152,157,159,163,165,167,169,173,175,178,181,184,186,190,192,195,197,201,204,207,209,209,210,222,223,226,227,232,231,238,237,244,243,247,247,251,252,257,258,259,262,267,268,269,270,269,278,281,282,285,286,289,290,293,294,294]
GROUP=[1,18,1,2,13,14,15,16,17,18,1,2,13,14,15,16,17,18]+[1,2]+list(range(3,19))+[1,2]+list(range(3,19))
streams={'element_symbols':''.join(SYM).upper(),'element_names_order':'HYDROGENHELIUMLITHIUMBERYLLIUMBORONCARBONNITROGENOXYGENFLUORINENEONSODIUMMAGNESIUMALUMINIUMSILICONPHOSPHORUSSULFURCHLORINEARGONPOTASSIUMCALCIUMSCANDIUMTITANIUMVANADIUMCHROMIUMMANGANESEIRONCOBALTNICKELCOPPERZINC'}
nums={'atomic_mass':MASS,'atomic_number':list(range(1,119)),'group':GROUP,
 'piano_freq_hz':[round(440*2**((k-49)/12)) for k in range(1,89)],'midi_C_major':[60,62,64,65,67,69,71,72]*20,
 'westminster':[64,60,62,55, 55,62,64,60, 64,62,60,55, 55,62,64,60]*10,'bach_motif':[70,69,72,71]*30,
 'chromatic_semitones':[k%12 for k in range(200)]}
def numscan(ct):
    best=0;hits=[]
    for nm,s in nums.items():
        ext=s*(200//len(s)+2)
        for an,A in (('ST',ST),('KA',KA)):
            c=np.array([A.index(ct[i]) for i in AN]);p=np.array([A.index(ANCH[i]) for i in AN])
            for mul in (1,-1):
                for rule in (0,1):
                    for off in range(len(s)):
                        k=np.array([(mul*ext[off+i])%26 for i in AN]);e=(p+k)%26 if rule==0 else (k-p)%26
                        m=int((e==c).sum());best=max(best,m)
                        if m==24: hits.append((nm,an,mul,rule,off))
    return hits,best
if __name__=='__main__':
    def parse(w):
        S=set(x.upper() for x in SYM);ok=[True]+[False]*len(w)
        for i in range(1,len(w)+1):
            ok[i]=any(ok[i-l] and w[i-l:i] in S for l in (1,2) if i-l>=0)
        return ok[-1]
    print('BERLINCLOCK spellable in element symbols:',parse('BERLINCLOCK'),'| EASTNORTHEAST:',parse('EASTNORTHEAST'),'| K4:',parse(CT0))
    h,b=numscan(CT0);print('numeric chemistry/music keys: hits',h,'best',b,'/24')
    l=list(CT0);random.Random(1).shuffle(l);print('null best',numscan(''.join(l))[1])
    big='|'.join(v+'|'+v[::-1] for v in streams.values())
    for an,A in (('KA',KA),('ST',ST)):
        arr,pos=prep(big*3,A);tot=0
        for W in range(1,9):
            Ms=np.arange(97,dtype=np.int16)[None,:] if W==1 else colperms(97,W)
            tot+=len(search(CT0,Ms,list(range(len(Ms))),arr,pos,A))
        print(an,'element-symbol / element-name streams as running key, direct + columnar 2-8: hits',tot)
