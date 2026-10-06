import cities as C
from cities import keyed,consistent_period,ST
eu={-1:'AZOREN KAPVERDE',0:'REYKJAVIK DUBLIN LONDON LISSABON MADEIRA BISSAU CASABLANCA CONAKRY DAKAR BAMAKO ACCRA',
1:'AMSTERDAM BERLIN BRUESSEL BRUSSEL BUDAPEST MADRID PARIS PRAG STOCKHOLM WARSCHAU OSLO KOPENHAGEN WIEN BERN PRESSBURG BELGRAD ROM TUNIS KINSHASA'}
top1='AMSTERDAM BERLIN BRUSSEL BUDAPEST MADRID PARIS PRAG STOCKHOLM WARSCHAU'.split()
bot1='OSLO KOPENHAGEN WIEN BERN PRESSBURG BELGRAD ROM TUNIS KINSHASA'.split()
words=[w for f in eu.values() for w in f.split()]+[''.join(top1),''.join(bot1),''.join(top1+bot1),''.join(eu[0].split()),''.join(eu[-1].split())]
# initials acrostics
words+=[''.join(w[0] for w in top1),''.join(w[0] for w in bot1)]
hits=[];n=0
for w in words:
    K=keyed(w)
    for name,A,B in (('Q1',K,ST),('Q2',ST,K),('Q3',K,K),('Q3rev',K[::-1],K[::-1])):
        for P in range(1,25):
            for rule in (0,1):
                n+=1;ok,red=consistent_period(A,B,P,rule)
                if ok and red>=1: hits.append((red,w,name,P,rule))
print('words',len(words),'tested',n,'hits',sorted(hits,reverse=True)[:10])
print('acrostics',words[-2:])
