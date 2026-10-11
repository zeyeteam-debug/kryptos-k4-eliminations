import time
from stack import CT0,ANCH
from sculptalpha import ALS
from biglist import WORDS
from layers3 import search,search_TTQ
from layers3b import search_QQT
sel=['PALIMPSEST','ABSCISSA','IQLUSION','UNDERGRUUND','DESPARATLY','BERLINCLOCK','EASTNORTHEAST','SHADOW','LUCID','DIGETAL','INTERPRETATU','LANGLEY','MAGNETIC','CANDLE','LAYERTWO','IDBYROWS']
t=time.time()
for w in sel:
    AL=ALS[w+'/ST'];h=search(CT0,ANCH,WORDS,AL);h2,h3=search_QQT(CT0,ANCH,WORDS,AL);h4=search_TTQ(CT0,ANCH,WORDS,AL)
    print(w,'TQQ',len(h),'QQT',len(h2),'QQQ',len(h3),'TTQ',len(h4),(h+h2+h3+h4)[:3],round(time.time()-t),flush=True)
