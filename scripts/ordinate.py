from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
def keyed(w):
    s=''
    for c in w+ST:
        if c not in s: s+=c
    return s
AN=sorted(ANCH)
words=['ORDINATE','ORDINATES','ORDINATA','YAXIS','COORDINATE','COORDINATES','ABSCISSAORDINATE','ORDINATEABSCISSA','APPLICATE','LATITUDE','LONGITUDE','ABSCISSA','PALIMPSEST','XAXIS','YCOORDINATE','NORTHING','EASTING']
best=[]
for w in words:
    al={'ST':ST,'KA':KA,'keyed(w)':keyed(w),'keyed(w+KRYPTOS)':keyed(w+'KRYPTOS'),'keyed(KRYPTOS+w)':keyed('KRYPTOS'+w)}
    for pn,PA in al.items():
        for cn,CA in al.items():
            for rule in range(3):
                for off in range(len(w)):
                    m=0
                    for i in AN:
                        k=CA.index(w[(i+off)%len(w)]);p=PA.index(ANCH[i])
                        m+=CA[[(p+k)%26,(k-p)%26,(p-k)%26][rule]]==CT0[i]
                    best.append((m,w,pn,cn,rule,off))
best.sort(reverse=True);print('combos',len(best),'best',best[:3])
