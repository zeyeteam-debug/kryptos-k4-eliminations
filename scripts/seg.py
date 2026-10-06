from stack import CT0,ANCH
import itertools
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
runs={'EASTNORTHEAST':range(21,34),'BERLINCLOCK':range(63,74)}
for nm,rg in runs.items():
    rg=list(rg);print('==',nm)
    for P in range(1,13):
        # any alphabets: same class -> same pt <=> same ct
        pairs=[(a,b) for a,b in itertools.combinations(rg,2) if a%P==b%P]
        strong=sum(ANCH[a]==ANCH[b] for a,b in pairs)
        anyok=all((ANCH[a]==ANCH[b])==(CT0[a]==CT0[b]) for a,b in pairs)
        vk=[]
        for an,A in (('ST',ST),('KA',KA)):
            for r in (0,1):
                K={i:((A.index(CT0[i])-A.index(ANCH[i])) if r==0 else (A.index(CT0[i])+A.index(ANCH[i])))%26 for i in rg}
                if all(K[a]==K[b] for a,b in pairs) and pairs: vk.append(an+('vig' if r==0 else 'beau'))
        print(f'P={P:2d} checks={len(pairs):2d} strong={strong} any-alphabet={"possible" if anyok else "NO"}  ST/KA vig/beau: {vk if vk else "NO"}')
