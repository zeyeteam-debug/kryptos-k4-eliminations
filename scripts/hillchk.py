import numpy as np,math,pickle,collections
from tauto import perms
from stack import CT0,ANCH
from hillT import run
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
Q,fl=pickle.load(open('quad.pkl','rb'))
def sc(t): return sum(Q.get(t[i:i+4],fl) for i in range(len(t)-3))/(len(t)-3)
M,Mi=perms(9)
h=run(CT0,[9],(2,),'K4')
inv_ok=0;results=[];dets=collections.Counter()
for hit in h:
    order,n,o,dname,x,k,an,W=hit
    m=(M if dname=='enc' else Mi)[x].astype(int)
    eqs=[(m[s],m[s+1],s) for s in range(o,96,2) if m[s] in ANCH and m[s+1] in ANCH]
    rows=[]
    for r in range(2):
        sols=[(a,b) for a in range(26) for b in range(26) if all((a*ST.index(ANCH[i])+b*ST.index(ANCH[j]))%26==ST.index(CT0[s+r]) for i,j,s in eqs)]
        rows.append(sols)
    for r0 in rows[0]:
        for r1 in rows[1]:
            det=(r0[0]*r1[1]-r0[1]*r1[0])%26;dets[math.gcd(det,26)]+=1
            if math.gcd(det,26)!=1: continue
            inv_ok+=1
            di=pow(det,-1,26);Inv=[[r1[1]*di%26,(-r0[1])*di%26],[(-r1[0])*di%26,r0[0]*di%26]]
            y=list(CT0)
            for s in range(o,96,2):
                c=[ST.index(CT0[s]),ST.index(CT0[s+1])]
                for r in range(2): y[s+r]=ST[(Inv[r][0]*c[0]+Inv[r][1]*c[1])%26]
            pt=['?']*97
            for j in range(97): pt[m[j]]=y[j]
            pt=''.join(pt);results.append((sc(pt),pt,hit))
print('hits',len(h),'matrix gcd(det,26) counts:',dets,'invertible:',inv_ok)
for r in sorted(results,reverse=True)[:5]: print(round(r[0],2),r[1])
