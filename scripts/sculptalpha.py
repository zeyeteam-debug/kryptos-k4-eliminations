# Keyed alphabets from words of the sculpture itself (K1-K3 plaintexts, keys, Morse panel, clues, misspellings)
import numpy as np,time,sys
from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
TXT="""BETWEEN SUBTLE SHADING AND THE ABSENCE OF LIGHT LIES THE NUANCE OF IQLUSION
IT WAS TOTALLY INVISIBLE HOWS THAT POSSIBLE THEY USED THE EARTHS MAGNETIC FIELD THE INFORMATION WAS GATHERED AND TRANSMITTED UNDERGRUUND
TO AN UNKNOWN LOCATION DOES LANGLEY KNOW ABOUT THIS THEY SHOULD ITS BURIED OUT THERE SOMEWHERE WHO KNOWS THE EXACT LOCATION ONLY WW THIS WAS HIS
LAST MESSAGE THIRTY EIGHT DEGREES FIFTY SEVEN MINUTES SIX POINT FIVE SECONDS NORTH SEVENTY SEVEN DEGREES EIGHT MINUTES FORTY FOUR SECONDS WEST ID BY ROWS
SLOWLY DESPARATLY SLOWLY THE REMAINS OF PASSAGE DEBRIS THAT ENCUMBERED THE LOWER PART OF THE DOORWAY WAS REMOVED WITH TREMBLING HANDS I MADE A TINY
BREACH IN THE UPPER LEFT HAND CORNER AND THEN WIDENING THE HOLE A LITTLE I INSERTED THE CANDLE AND PEERED IN THE HOT AIR ESCAPING FROM THE CHAMBER
CAUSED THE FLAME TO FLICKER BUT PRESENTLY DETAILS OF THE ROOM WITHIN EMERGED FROM THE MIST CAN YOU SEE ANYTHING
VIRTUALLY INVISIBLE DIGETAL INTERPRETATU SHADOW FORCES LUCID MEMORY T IS YOUR POSITION SOS RQ
PALIMPSEST ABSCISSA KRYPTOS BERLIN CLOCK BERLINCLOCK EAST NORTHEAST EASTNORTHEAST LAYERTWO IDBYROWS WHATSTHEPOINT DELIVERINGAMESSAGE
ITSBURIEDOUTTHERESOMEWHERE EGYPT NOVEMBER WALL"""
def keyed(w,base=ST):
    s=''
    for c in w+base:
        if c not in s: s+=c
    return s
words=sorted({w for w in TXT.split() if len(w)>=3})
ALS={};seen=set()
for w in words:
    for base,bn in ((ST,'ST'),(ST[::-1],'STrev')):
        a=keyed(w,base)
        if a in seen: continue
        seen.add(a);ALS[w+'/'+bn]=a
if __name__=='__main__':
    t=time.time();print('alphabets',len(ALS),flush=True)
    import dictsearch as D
    # (1) Q(370k dictionary) -> none / 752 spirals ; and Q -> columnar 2..7 (all orders), both orders
    sp=[]
    for Wd in range(2,49):
        for corner in ('TL','TR','BL','BR'):
            for cw in (0,1):
                m=D.perm_from_read(D.spiral_order(97,Wd,corner,cw),97);sp+=[((Wd,corner,cw,'r'),m),((Wd,corner,cw,'w'),np.argsort(m))]
    tot=0
    for k,(nm,AL) in enumerate(ALS.items()):
        TB=D.table(D.WORDS,ANCH,AL)
        h=D.run(CT0,ANCH,sp,[],AL,TB)+D.run(CT0,ANCH,[('none',np.arange(97))],[2,3,4,5,6,7],AL,TB)
        tot+=len(h)
        if h: print('HIT',nm,h[:5],flush=True)
        if k%20==0: print(k,nm,round(time.time()-t),flush=True)
    print('(1) total hits',tot,flush=True)
