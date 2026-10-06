import numpy as np,glob,os,sys,random,time
from stack import CT0,ANCH
from gen import colperms
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
AN=np.array(sorted(ANCH))
def load_texts():
    T={}
    base='texts/nltk_data/packages/corpora'
    for f in sorted(glob.glob(base+'/gutenberg/*.txt')+glob.glob(base+'/inaugural/*.txt')+glob.glob(base+'/state_union/*.txt')):
        s=''.join(ch for ch in open(f,errors='ignore').read().upper() if 'A'<=ch<='Z')
        T[os.path.basename(f)]=s
    exec(open('board.py').read().split('# K4 cell')[0],globals())
    T['K1K2K3']=K1P+K2P+K3P
    big='|'.join(T.values())        # separators break matches across files
    return T,big
def prep(big,A):
    lut=np.full(128,255,dtype=np.uint8)
    for i,c in enumerate(A): lut[ord(c)]=i
    arr=lut[np.frombuffer(big.encode(),dtype=np.uint8)]
    pos=[np.nonzero(arr==v)[0] for v in range(26)]
    return arr,pos
def search(ct,Ms,labs,arr,pos,A,label=''):
    C=np.array([A.index(x) for x in ct]);P=np.array([A.index(ANCH[i]) for i in AN]);L=len(arr);hits=[]
    for x in range(len(Ms)):
        q=Ms[x][AN]
        for rule in (0,1,2):
            d=((C[q]-P) if rule==0 else (C[q]+P) if rule==1 else (P-C[q]))%26
            for kind,idx in (('ct',q),('pt',AN)):
                o=np.argsort(idx);ii=idx[o];dd=d[o]                  # anchor positions relative to key start
                cand=pos[dd[0]]-ii[0]
                cand=cand[(cand>=0)&(cand+ii.max()<L)]
                for j in range(1,len(ii)):
                    if len(cand)==0: break
                    cand=cand[arr[cand+ii[j]]==dd[j]]
                for off in cand: hits.append((labs[x],rule,kind,int(off)))
    return hits
if __name__=='__main__':
    Ws=[int(w) for w in sys.argv[1].split(',')]
    T,big=load_texts();print('text letters:',len(big),'files',len(T),flush=True)
    for an,A in (('KA',KA),('ST',ST)):
        arr,pos=prep(big,A)
        # control: Moby Dick key at offset 123456, columnar width 6, key after transposition
        if an=='KA':
            src='THEYUSEDTHEEARTHSMAGNETICFIELDXTHEINFORMATIONWASGATHEREDANDTRANSMITTEDUNDERGROUNDTOANUNKNOWNLOCATIONXDOES'
            pt=list(src[:97]);pt[21:34]='EASTNORTHEAST';pt[63:74]='BERLINCLOCK';pt=''.join(pt)
            M6=colperms(97,6);m=M6[100];y=['?']*97
            for i in range(97): y[m[i]]=pt[i]
            off=123456;ctrl=''.join(KA[(KA.index(c)+int(arr[off+j]))%26] for j,c in enumerate(y))
            t=time.time();h=search(ctrl,M6[100:101],['ctrl'],arr,pos,KA);print('control:',h,round(time.time()-t,2),'s/perm',flush=True)
        for W in Ws:
            Ms=colperms(97,W);t=time.time()
            h=search(CT0,Ms,[(W,i) for i in range(len(Ms))],arr,pos,A)
            print(an,'W',W,'perms',len(Ms),'hits',len(h),h[:5],round(time.time()-t),'s',flush=True)
