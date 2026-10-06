import itertools
from stack import solve_mod,ANCH,CT0
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
def keyed(w):
    s=''
    for c in w+ST:
        if c not in s: s+=c
    return s
faces={5:'JEKATERINBURG BISCHKEK DUSCHANBE NEWDELHI KARACHI COLOMBO',6:'OMSK ALMATY TASCHKENT NOWOSIBIRSK RANGUN DHAKA',7:'KRASNOJARSK HANOI BANGKOK PHNOMPENH JAKARTA',8:'IRKUTSK ULANBATOR PEKING SHANGHAI MANILA PERTH HONGKONG KUALALUMPUR SINGAPUR',9:'JAKUTSK PJONGJANG TOKYO SEOUL',10:'CHABAROWSK WLADIWOSTOK SYDNEY CANBERRA MELBOURNE',11:'MAGADAN SACHALIN',12:'KAMTSCHATKA DATUMSGRENZE WELLINGTON'}
extra='BERLIN MOSKAU MOSCOW WELTZEITUHR URANIA ALEXANDERPLATZ KRYPTOS'.split()
words=[w for f in faces.values() for w in f.split()]+extra
# also each full face concatenated, and the sequence of faces
words+=[''.join(f.split()) for f in faces.values()]
idx=sorted(ANCH)
def consistent_period(A,B,P,rule):
    # Quagmire: ct = B[(A^-1? )]: use model ct_index_in_B = pt_index_in_A + k[i%P]  (rule0) or k - pt (rule1)
    rows=[];rhs=[]
    for i in idx:
        row=[0]*P;row[i%P]=1
        c=B.index(CT0[i]);p=A.index(ANCH[i])
        rows.append(row);rhs.append((c-p)%26 if rule==0 else (c+p)%26)
    ok2,_=solve_mod(rows,rhs,2);ok13,r=solve_mod(rows,rhs,13)
    return ok2 and ok13, len(idx)-r
hits=[];n=0
for w in words:
    K=keyed(w)
    for name,A,B in (('Q1',K,ST),('Q2',ST,K),('Q3',K,K),('Q3rev',K[::-1],K[::-1])):
        for P in range(1,25):
            for rule in (0,1):
                n+=1;ok,red=consistent_period(A,B,P,rule)
                if ok and red>=1: hits.append((red,w,name,P,rule))
print('tested',n,'consistent:',sorted(hits,reverse=True)[:10])
