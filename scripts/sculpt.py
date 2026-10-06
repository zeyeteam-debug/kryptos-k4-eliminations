import math
from reportlab.platypus import Flowable
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
pdfmetrics.registerFont(TTFont('DVMB','/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf'))
KA='KRYPTOSABCDEFGHIJLMNQUVWXZ';ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
HDR=' ABCDEFGHIJKLMNOPQRSTUVWXYZABCD'
TAB=[HDR]+[ST[r]+''.join(KA[(r+j)%26] for j in range(30))+('L' if ST[r]=='N' else '') for r in range(26)]+[HDR]
class Sculpture(Flowable):
    def __init__(self,lines,width):
        super().__init__();self.lines=lines;self.width=width
        self.cols=33+3+32;self.adv=width/(self.cols+3);self.fs=self.adv/0.6;self.lh=self.fs*1.12
        self.height=len(lines)*self.lh+self.lh*1.6
    def draw(self):
        c=self.canv;W,H=self.width,self.height
        # copper sheet with S-curve shading
        n=160
        for k in range(n):
            x=W*k/n;t=k/n
            sh=0.5+0.5*math.cos(2*math.pi*(t*1.5+0.1))     # light/dark bands of an S-shaped screen
            r=0.36+0.20*sh;g=0.25+0.15*sh;b=0.15+0.08*sh
            c.setFillColorRGB(r,g,b);c.rect(x,0,W/n+0.6,H,stroke=0,fill=1)
        c.setStrokeColorRGB(0.22,0.15,0.09);c.setLineWidth(0.6)
        top=H-self.lh*1.05
        seam_y=top-14*self.lh+self.lh*0.25
        c.line(0,seam_y,W,seam_y)                       # top/bottom plates
        xL=self.adv;xR=self.adv*(1+33+3)
        c.line(xR-self.adv*1.5,0,xR-self.adv*1.5,H)    # left (cipher) / right (tableau)
        c.setFillColorRGB(0.08,0.05,0.03)
        def put(x0,y,s,raised=()):
            for i,ch in enumerate(s):
                if ch==' ': continue
                dy=self.fs*0.35 if i in raised else 0
                c.setFont('DVMB',self.fs*(0.8 if i in raised else 1));c.drawString(x0+i*self.adv,y+dy,ch)
        for L,line in enumerate(self.lines):
            y=top-L*self.lh-(self.lh*0.5 if L>=14 else 0)
            put(xL,y,line,raised=(3,4,6) if L==14 else ())
        for L,line in enumerate(TAB):
            y=top-L*self.lh-(self.lh*0.5 if L>=14 else 0)
            put(xR,y,line)
