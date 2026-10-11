from stack import CT0,ANCH
ST='ABCDEFGHIJKLMNOPQRSTUVWXYZ';KA='KRYPTOSABCDEFGHIJLMNQUVWXZ'
AN=sorted(ANCH);res=[];n=0
for an,A in (('ST',ST),('KA',KA)):
  for base in (0,1):
    for pb in ('id','swap21'):
      def code(ch):
        v=A.index(ch)+base
        if pb=='swap21':
          b=[(v>>i)&1 for i in range(5)];b[1],b[2]=b[2],b[1];v=sum(x<<i for i,x in enumerate(b))
        return v
      k={i:code(CT0[i])^code(ANCH[i]) for i in AN}
      for P in range(1,25):
        n+=1;d={};ok=True
        for i in AN:
          if d.setdefault(i%P,k[i])!=k[i]: ok=False;break
        red=24-len({i%P for i in AN})
        if ok and red>=1: res.append((an,base,pb,P,red))
print('5-bit XOR periodic keys tested',n,'consistent:',res)
