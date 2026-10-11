# 5-bit / XOR (Vernam) family: codes A0, A1 (ST/KA order) and ITA2 (Baudot-Murray teleprinter code)
import numpy as np,random
from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
ITA2=dict(A='11000',B='10011',C='01110',D='10010',E='10000',F='10110',G='01011',H='00101',I='01100',J='11010',K='11110',L='01001',M='00111',
 N='00110',O='00011',P='01101',Q='11101',R='01010',S='10100',T='00001',U='11100',V='01111',W='11001',X='10111',Y='10101',Z='10001')
AN=sorted(ANCH)
def encs():
    E={}
    for an,A in (('ST',ST),('KA',KA)):
        for b in (0,1): E[f'{an}{b}']={c:A.index(c)+b for c in ST}
    E['ITA2']={c:int(v,2) for c,v in ITA2.items()}
    E['ITA2rev']={c:int(v[::-1],2) for c,v in ITA2.items()}
    return E
E=encs()
def keyvals(ct,enc):
    return [enc[ct[i]]^enc[ANCH[i]] for i in AN]
def periodic(ct):
    out=[]
    for nm,enc in E.items():
        k=dict(zip(AN,keyvals(ct,enc)))
        for P in range(1,25):
            d={};ok=all(d.setdefault(i%P,k[i])==k[i] for i in AN);red=24-len({i%P for i in AN})
            if ok and red>=1: out.append((nm,P,red))
    return out
def autokey(ct):
    out=[]
    for nm,enc in E.items():
        for L in range(1,41):
            for kind in ('pt','ct'):
                chk=0;ok=True
                for i in AN:
                    src=ANCH.get(i-L) if kind=='pt' else (ct[i-L] if i-L>=0 else None)
                    if src is None: continue
                    chk+=1
                    if enc[ct[i]]!=enc[ANCH[i]]^enc[src]: ok=False;break
                if ok and chk>=4: out.append((nm,L,kind,chk))
    return out
def running(ct,T):
    out=[];g=42
    for nm,enc in E.items():
        inv={v:c for c,v in enc.items()};kv=keyvals(ct,enc)
        if any(v not in inv for v in kv): out.append((nm,'impossible: XOR gives a non-letter code'));continue
        key=''.join(inv[v] for v in kv);f1,f2=key[:13],key[13:]
        for tn,txt in T.items():
            for rev in (0,1):
                s=txt[::-1] if rev else txt;st=s.find(f1)
                while st!=-1:
                    if s[st+g:st+g+11]==f2: out.append((nm,tn,rev,st))
                    st=s.find(f1,st+1)
    return out
def bm(s):
    # Berlekamp-Massey over GF(2); returns linear complexity
    n=len(s);c=[1]+[0]*n;b=[1]+[0]*n;L=0;m=-1
    for N_ in range(n):
        d=s[N_]
        for i in range(1,L+1): d^=c[i]&s[N_-i]
        if d:
            t=c[:];
            for i in range(n-N_+m):
                if N_-m+i<=n: c[N_-m+i]^=b[i]
            if L<=N_//2: L=N_+1-L;m=N_;b=t
    return L
def lfsr(ct):
    out=[]
    for nm,enc in E.items():
        for order in ('msb','lsb'):
            for blk in ((21,34),(63,74)):
                bits=[]
                for i in range(*blk):
                    v=enc[ct[i]]^enc[ANCH[i]];bb=[(v>>(4-j))&1 for j in range(5)]
                    bits+=bb if order=='msb' else bb[::-1]
                out.append((nm,order,blk,len(bits),bm(bits)))
    return out
if __name__=='__main__':
    SH=list(CT0);random.Random(7).shuffle(SH);SH=''.join(SH)
    # controls
    pt=list('THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'[:97])
    for i,c in ANCH.items(): pt[i]=c
    enc=E['ITA2'];inv={v:c for c,v in enc.items()}
    # ct-autokey control in ITA2 needs XOR results to be letter codes; plant with LFSR keystream bits instead
    st=[1,0,0,1,1,0,1,0,1,1,1,0,0,0,1,0,1];taps=(17,14)
    bits=[]
    s=st[:]
    for _ in range(97*5):
        nb=s[-taps[0]]^s[-taps[1]];bits.append(s[-1]);s.append(nb)
    print('LFSR control: LC of planted keystream block =',bm(bits[21*5:34*5]),'(degree 17)')
    for nm,ct in (('K4',CT0),('null',SH)):
        print(nm,'periodic XOR:',periodic(ct))
        print(nm,'autokey XOR:',autokey(ct))
        L=lfsr(ct);print(nm,'LFSR linear complexity (block bits, LC): min',min(x[4] for x in L),'| all',[(x[0],x[1],x[2][0],x[3],x[4]) for x in L][:6])
    from allalpha import texts
    T=texts()
    for nm,ct in (('K4',CT0),('null',SH)): print(nm,'running-key XOR:',running(ct,T)[:6])
