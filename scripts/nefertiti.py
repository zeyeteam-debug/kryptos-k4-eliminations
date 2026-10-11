import numpy as np,random,time
import allalpha as AA
from stack import CT0,ANCH
NW="""NEFERTITI NOFRETETE NOFERTITI NEFERNEFERUATEN NEFERNEFERUATENNEFERTITI AKHENATEN AKHENATON ECHNATON AMENHOTEP AMENOPHIS ATEN ATON
AMARNA TELLELAMARNA AKHETATEN THUTMOSE THUTMOSIS BORCHARDT LUDWIGBORCHARDT JAMESSIMON SIMON BUST QUEEN QUEENNEFERTITI THEBEAUTIFULONEHASCOME
CHARLOTTENBURG STUELERBAU AEGYPTISCHESMUSEUM EGYPTIANMUSEUM NEUESMUSEUM MUSEUMSINSEL TUTANKHAMUN TUTANKHAMEN TUTANKHATEN SMENKHKARE
MERITATEN ANKHESENAMUN KIYA HOWARDCARTER CARTER CARNARVON VALLEYOFTHEKINGS KV LUXOR THEBES KARNAK GIZA CAIRO SAKKARA PYRAMID SPHINX
NILE HIEROGLYPH ROSETTA CARTOUCHE SCARAB ANKH EYEOFHORUS HORUS RA OSIRIS ISIS ANUBIS HATHOR PHARAOH DYNASTY EIGHTEENTH ONEEYE""".split()
NW=sorted(set(NW))
def keyed(w,base):
    return ''.join(dict.fromkeys(w+base))
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
ALS={}
for w in NW:
    for b,bn in ((ST,'ST'),(ST[::-1],'STrev'),(KA,'KA')):
        a=keyed(w,b)
        if a not in ALS.values(): ALS[w+'/'+bn]=a
if __name__=='__main__':
    t=time.time();SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH);print('words',len(NW),'alphabets',len(ALS),flush=True)
    AA.ALS=ALS;M=AA.linear_models()
    for nm,ct in (('K4',CT0),('null',SH)):
        labs,K=AA.Kmat(ct)
        print(nm,'linear',len(AA.linear_test(K,M)),'safe',len(AA.safe_test(K)),'masks',len(AA.mask_test(K)[0]),'hill',len(AA.hill_test(ct)),'autokey',len(AA.autokey_test(ct)),round(time.time()-t),flush=True)
    T=AA.texts()
    for nm,ct in (('K4',CT0),('null',SH)): print(nm,'running key',len(AA.rk_test(ct,T)),flush=True)
    import dictsearch as D
    sp=[('none',np.arange(97))]
    for Wd in range(2,49):
        for corner in ('TL','TR','BL','BR'):
            for cw in (0,1):
                m=D.perm_from_read(D.spiral_order(97,Wd,corner,cw),97);sp+=[((Wd,corner,cw,'r'),m),((Wd,corner,cw,'w'),np.argsort(m))]
    tot=0
    # (a) Nefertiti words as Q-keys with ST/KA base, spirals + columnar 2-7
    for an,A in (('ST',ST),('KA',KA)):
        TB=D.table(NW,ANCH,A);h=D.run(CT0,ANCH,sp,[2,3,4,5,6,7],A,TB);tot+=len(h);print('(a) Q(Nefertiti words)->spiral->columnar',an,len(h),flush=True)
    # (b) Nefertiti-keyed alphabets with the full 370k dictionary as Q key
    for nm,A in ALS.items():
        TB=D.table(D.WORDS,ANCH,A);h=D.run(CT0,ANCH,sp,[],A,TB)+D.run(CT0,ANCH,[('none',np.arange(97))],[2,3,4,5,6,7],A,TB);tot+=len(h)
        if h: print('HIT',nm,h[:3],flush=True)
    print('(b) dictionary Q with Nefertiti alphabets: total',tot,round(time.time()-t),flush=True)
    from biglist import WORDS
    from layers3 import search,search_TTQ
    from layers3b import search_QQT
    W=WORDS+NW
    for an,A in (('ST',ST),('KA',KA),('NEFERTITI',ALS['NEFERTITI/ST']),('AKHENATEN',ALS['AKHENATEN/ST'])):
        h=search(CT0,ANCH,W,A);h2,h3=search_QQT(CT0,ANCH,W,A);h4=search_TTQ(CT0,ANCH,W,A)
        print('(c) layers',an,len(h),len(h2),len(h3),len(h4),round(time.time()-t),flush=True)
