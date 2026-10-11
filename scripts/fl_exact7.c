// Exact DFS over 10x10 Fleissner grilles (4^25) + periodic substitution consistency on the K4 cribs.
// stdin: line1 ct(97) ; line2 dir(+1 cw / -1 ccw) blankmode(0=end,1=start) narb arbp...;
//        line3: 100 ints write order (cell ids r*10+c) ; line4: 100 ints read order
// Two-ended processing: step s handles write-order cell s (front) and 99-s (back) -- same orbit (180 deg).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static const char *KA="KRYPTOSABCDEFGHIJLMNQUVWXZ",*ST="ABCDEFGHIJKLMNOPQRSTUVWXYZ";
static char CT[128];int dir,bmode,narb,arbp[16];
int W[100],RD[100],wpos[100],rpos[100];
int orb[100],oj[100],ocell[25][4];
int lab[100];int branched[25];
int crib[100];// k4 index -> letter or -1
int cidx[2][128];
// models: 48 vig: m = a*24 + t*12 + (p-1)
int val[48][12];
int amapP[16][12][26],amapC[16][12][26];int vmask[48][12];int pmask[16][12][26];
uint64_t alive;
int *ulog_p[1<<16];int ulog_v[1<<16];int ulen;
int cntF[4],cntB[4];
int blanks[3],nbl;int pend[100][2],npend;
int cur_s;long long nodes,leaves,dn[51],dck[51];double frac,wt=1.0;int nbr;double t0;
#include <time.h>
static double now(){struct timespec ts;clock_gettime(CLOCK_MONOTONIC,&ts);return ts.tv_sec+1e-9*ts.tv_nsec;}
int labsol[100];
FILE *out;
static inline void setv(int *p,int v){ulog_p[ulen]=p;ulog_v[ulen]=*p;ulen++;*p=v;}
static void undo_to(int L){while(ulen>L){ulen--;*ulog_p[ulen]=ulog_v[ulen];}}
// add crib constraint: k4 index i, ciphertext index c. returns 0 if all models dead
unsigned char modtab[48][100],dtab[48][26][26];
static int addc(int i,int c){
  int C=CT[c]-'A',P=crib[i]-'A';
  uint64_t a=alive&((1ULL<<48)-1);
  while(a){int m=__builtin_ctzll(a);a&=a-1;
    int d=dtab[m][C][P];int *v=&val[m][modtab[m][i]];
    if(*v<0){ if(vmask[m][modtab[m][i]]>>d&1) setv(v,d); else alive&=~(1ULL<<m);} else if(*v!=d) alive&=~(1ULL<<m);
  }
  for(int k=0;k<narb;k++) if(alive>>(48+k)&1){
    int p=arbp[k],cl=i%p,pp=P,cc=C;
    int *x=&amapP[k][cl][pp],*y=&amapC[k][cl][cc];
    if(!(pmask[k][cl][pp]>>cc&1)) {alive&=~(1ULL<<(48+k));continue;}
    if(*x<0 && *y<0){setv(x,cc);setv(y,pp);}
    else if(*x!=cc || *y!=pp) alive&=~(1ULL<<(48+k));
  }
  return alive!=0;
}

// relaxed constraint: crib i lies at one of the ciphertext indices cs[0..n-1]
static int relax(int i,int *cs,int n){
  int P=crib[i]-'A';
  uint64_t a=alive&((1ULL<<48)-1);
  while(a){int m=__builtin_ctzll(a);a&=a-1;
    int bits=0;for(int j=0;j<n;j++) bits|=1<<dtab[m][CT[cs[j]]-'A'][P];
    int cl=modtab[m][i];int v=val[m][cl];
    if(v>=0){ if(!(bits>>v&1)) alive&=~(1ULL<<m);} 
    else { int nm=vmask[m][cl]&bits; if(!nm) alive&=~(1ULL<<m); else if(nm!=vmask[m][cl]) setv(&vmask[m][cl],nm);}
  }
  for(int k=0;k<narb;k++) if(alive>>(48+k)&1){
    int p=arbp[k],cl=i%p;int bits=0,any=0;
    for(int j=0;j<n;j++){int cc=CT[cs[j]]-'A';bits|=1<<cc; if(amapC[k][cl][cc]<0||amapC[k][cl][cc]==P) any=1;}
    int x=amapP[k][cl][P];
    if(x>=0){ if(!(bits>>x&1)) alive&=~(1ULL<<(48+k)); continue;}
    int nm=pmask[k][cl][P]&bits; if(!nm||!any) alive&=~(1ULL<<(48+k)); else if(nm!=pmask[k][cl][P]) setv(&pmask[k][cl][P],nm);
  }
  return alive!=0;
}
static int ctindex(int x){int r=rpos[x],n=0;for(int k=0;k<3;k++) if(rpos[blanks[k]]<r) n++;return r-n;}
// handle a cell whose plaintext index is now known; return 0 to prune
int tied,preset,nloc;long long fams;double famcount;int req[100],forb[100];
static int handle(int x,int pt,int bb){
  int isb = bmode==0 ? pt>=97 : pt<3;
  if(isb){ if(preset) return 1; blanks[nbl++]=x;
    if(nbl==3){ for(int k=0;k<npend;k++) if(!addc(pend[k][0],ctindex(pend[k][1]))) return 0; }
    return 1;}
  int i = bmode==0 ? pt : pt-3;
  if(crib[i]<0) return 1;
  nloc++;
  if(tied){dck[cur_s]++;return addc(i,rpos[x]-bb);}
  if(nbl<3){pend[npend][0]=i;pend[npend][1]=x;npend++;return 1;}
  dck[cur_s]++;return addc(i,ctindex(x));
}
static void dfs(int s);
int lastg=-1,lam;
static int isblank(int x){for(int k=0;k<nbl;k++) if(blanks[k]==x) return 1;return 0;}
// shift constraints for cribs placed before all blanks are known: exactly k blanks precede read position c
int shc[64],shk[64],nsh;
static int known_before(int c){int n=0;for(int k=0;k<nbl;k++) if(rpos[blanks[k]]<c) n++;return n;}
static int shifts_ok(void){int rem=3-nbl;
  for(int j=0;j<nsh;j++){int kb=known_before(shc[j]);if(kb>shk[j]||kb+rem<shk[j]) return 0;}
  return 1;}
static int future_lam(int s){int n=0;
  for(int t=(lastg>s?lastg:s)+1;t<50;t++){int y=W[t];int l=lab[y];
    int can = l>=0 ? l==lam : !(forb[y]>>lam&1) && (req[y]<0||req[y]==lam); n+=can;} return n;}
static void place(int s,int *ci,int *cx,int nc,int reg){
  int ok=1;
  if(nbl<3){
    if(future_lam(s)<3-nbl) ok=0;   // top half would hold <3 label-lam cells: covered by the rare-case pre-branch
    for(int k=0;k<nc&&ok;k++){pend[npend][0]=ci[k];pend[npend][1]=cx[k];npend++;
      int c=rpos[cx[k]],nbefore=0,nafter=0;
      for(int t=(lastg>s?lastg:s)+1;t<50;t++){int y=W[t];int l=lab[y];
        int can = l>=0 ? l==lam : !(forb[y]>>lam&1) && (req[y]<0||req[y]==lam);
        if(!can) continue; int z= bmode==0 ? W[99-t] : y; if(rpos[z]<c) nbefore++; else nafter++;}
      int rem=3-nbl,kb=known_before(c);
      int kmin=kb+(rem-nafter>0?rem-nafter:0),kmax=kb+(rem<nbefore?rem:nbefore);
      int cs[4],n=0;for(int kk=kmin;kk<=kmax;kk++) cs[n++]=c-kk;
      if(n==0) ok=0; else {dck[s]++;ok=relax(ci[k],cs,n);} }
  } else {
    if(reg) for(int k=0;k<npend&&ok;k++){dck[s]++;ok=addc(pend[k][0],ctindex(pend[k][1]));}
    for(int k=0;k<nc&&ok;k++){dck[s]++;ok=addc(ci[k],ctindex(cx[k]));}
  }
  if(ok) dfs(s+1); else frac+=wt;
}
static void process(int s){
  int f=W[s],b=W[99-s];cur_s=s;
  int L0=ulen;uint64_t al=alive;int nb=nbl,np=npend,nl=nloc,lg=lastg;
  int bbF = bmode==0 ? (cntF[3]>22?cntF[3]-22:0) : (cntF[0]<3?cntF[0]:3);
  int bbB = bmode==0 ? 3-(cntB[3]<3?cntB[3]:3) : 3-(cntB[0]>22?cntB[0]-22:0);
  int tf=lab[f];int ptf=25*tf+cntF[tf];cntF[tf]++;
  int tb=lab[b];int ptb=25*tb+24-cntB[tb];cntB[tb]++;
  if(tied){
    if(handle(f,ptf,bbF) && handle(b,ptb,bbB)) dfs(s+1); else frac+=wt;
  } else {
    int xs[2]={f,b},ps[2]={ptf,ptb};
    int ci[2],cx[2],nc=0,reg=0;
    for(int k=0;k<2;k++){int pt=ps[k];int isb = bmode==0 ? pt>=97 : pt<3;
      if(isb){ if(!isblank(xs[k])){ blanks[nbl++]=xs[k];reg=1; } if(s>lastg) lastg=s; continue;}
      int i = bmode==0 ? pt : pt-3; if(crib[i]>=0){ci[nc]=i;cx[nc]=xs[k];nc++;nloc++;}}
    place(s,ci,cx,nc,reg);
  }
  cntF[tf]--;cntB[tb]--;undo_to(L0);alive=al;nbl=nb;npend=np;nloc=nl;lastg=lg;
}
static void dfs(int s){
  nodes++;dn[s]++;
  if((nodes&((1<<28)-1))==0){fprintf(stderr,"prog nodes=%lld frac=%.6g t=%.0f\n",nodes,frac,now()-t0);for(int q=0;q<51;q++)fprintf(stderr,"%d:%lld/%lld ",q,dn[q],dck[q]);fprintf(stderr,"\n");}
  if(nloc==24 && nbl==3 && s<50){
    // all cribs and blanks located and consistent: every completion of the free orbits survives
    int fr=0;for(int o=0;o<25;o++) if(!branched[o]) fr++;
    fams++;famcount+=(double)(1ULL<<(2*fr));frac+=wt;
    fprintf(out,"FAM alive=%llx free=%d lab=",(unsigned long long)alive,fr);
    for(int x=0;x<100;x++) fprintf(out,"%c",lab[x]<0?'.':'0'+lab[x]);
    fprintf(out,"\n");fflush(out);return;}
  if(s==50){ leaves++;frac+=wt;famcount+=1;
    // record grille: H0 = cells with label 0
    fprintf(out,"SOL alive=%llx H0=",(unsigned long long)alive);
    for(int x=0;x<100;x++) if(lab[x]==0) fprintf(out,"%d,",x);
    fprintf(out,"\n");fflush(out);return;}
  int f=W[s],o=orb[f];
  if(branched[o]){process(s);return;}
  branched[o]=1;wt/=4;
  for(int L=0;L<4;L++){
    int ok=1;
    for(int k=0;k<4;k++){int d=(k-oj[f]+4)%4;int y=ocell[o][k],l=((L+dir*d)%4+4)%4;lab[y]=l;
      if((req[y]>=0&&req[y]!=l)||(forb[y]>>l&1)) ok=0;}
    if(ok) process(s);
  }
  for(int k=0;k<4;k++) lab[ocell[o][k]]=-1;
  branched[o]=0;wt*=4;
}
int main(int argc,char**argv){
  out=stdout;
  if(scanf("%127s",CT)!=1) return 1;
  (void)!scanf("%d %d %d",&dir,&bmode,&narb);for(int k=0;k<narb;k++) scanf("%d",&arbp[k]);
  for(int k=0;k<100;k++) scanf("%d",&W[k]);
  for(int k=0;k<100;k++) scanf("%d",&RD[k]);
  for(int k=0;k<100;k++){wpos[W[k]]=k;rpos[RD[k]]=k;}
  tied=1;for(int k=0;k<100;k++) if(W[k]!=RD[k]) tied=0;
  // orbits: o_{j+1}=R(o_j), R=(r,c)->(c,9-r) clockwise
  int seen[100]={0},no=0;
  for(int x=0;x<100;x++) if(!seen[x]){int y=x;for(int j=0;j<4;j++){ocell[no][j]=y;orb[y]=no;oj[y]=j;seen[y]=1;int r=y/10,c=y%10;y=c*10+(9-r);}no++;}
  for(int i=0;i<100;i++) crib[i]=-1;
  const char*w1="EASTNORTHEAST",*w2="BERLINCLOCK";
  for(int j=0;w1[j];j++) crib[21+j]=w1[j];
  for(int j=0;w2[j];j++) crib[63+j]=w2[j];
  for(int k=0;k<26;k++){cidx[0][(int)ST[k]]=k;cidx[1][(int)KA[k]]=k;}
  for(int m=0;m<48;m++){int a=m/24,t=(m/12)%2,p=m%12+1;for(int i=0;i<100;i++)modtab[m][i]=i%p;
    for(int C=0;C<26;C++)for(int P=0;P<26;P++){int cc=cidx[a]['A'+C],pp=cidx[a]['A'+P];dtab[m][C][P]= t==0?(cc-pp+26)%26:(cc+pp)%26;}}
  memset(val,-1,sizeof val);for(int m=0;m<48;m++)for(int c=0;c<12;c++)vmask[m][c]=(1<<26)-1;for(int k=0;k<16;k++)for(int c=0;c<12;c++)for(int q=0;q<26;q++)pmask[k][c][q]=(1<<26)-1;memset(amapP,-1,sizeof amapP);memset(amapC,-1,sizeof amapC);
  alive=(1ULL<<(48+narb))-1;
  if(argc>1) alive&=strtoull(argv[1],0,16);
  t0=now();
  for(int x=0;x<100;x++){lab[x]=-1;req[x]=-1;forb[x]=0;}
  if(tied){ dfs(0); }
  else {
    lam = bmode==0 ? 1 : 0;
    preset=0;nbl=0;dfs(0);
    fprintf(stderr,"main sweep done nodes=%lld t=%.0f\n",nodes,now()-t0);
    // rare case: top half holds only m<3 label-lam cells (Q); the other 3-m blanks are the last 3-m
    // label-mu cells U of the top half (end: mu=3, blank=U itself; start: mu=2, blank=180-deg partner of U)
    int mu = bmode==0 ? 3 : 2;
    for(int m=0;m<3;m++){
      int nq = m==0?1:(m==1?50:2500), nu = (3-m)==3?125000:((3-m)==2?2500:50);
      for(int zq=0;zq<nq;zq++){
        int Q[2]; if(m==1)Q[0]=zq; if(m==2){Q[0]=zq/50;Q[1]=zq%50;if(Q[0]>=Q[1]) continue;}
        for(int zu=0;zu<nu;zu++){
          int U[3],k=3-m,okk=1;
          if(k==1)U[0]=zu; if(k==2){U[0]=zu/50;U[1]=zu%50;if(U[0]>=U[1]) continue;}
          if(k==3){U[0]=zu/2500;U[1]=(zu/50)%50;U[2]=zu%50;if(U[0]>=U[1]||U[1]>=U[2]) continue;}
          for(int j=0;j<m;j++)for(int i=0;i<k;i++) if(Q[j]==U[i]) okk=0;
          if(!okk) continue;
          for(int x=0;x<100;x++){req[x]=-1;forb[x]=0;}
          for(int t=0;t<50;t++) forb[W[t]]|=1<<lam;
          for(int t=U[0]+1;t<50;t++) forb[W[t]]|=1<<mu;
          for(int j=0;j<m;j++){forb[W[Q[j]]]&=~(1<<lam);req[W[Q[j]]]=lam;}
          for(int i=0;i<k;i++){forb[W[U[i]]]&=~(1<<mu);req[W[U[i]]]=mu;}
          int nb_=0;
          for(int j=0;j<m;j++) blanks[nb_++]= bmode==0 ? W[99-Q[j]] : W[Q[j]];
          for(int i=0;i<k;i++) blanks[nb_++]= bmode==0 ? W[U[i]] : W[99-U[i]];
          preset=1;nbl=3;npend=0;dfs(0);
        }
      }
      fprintf(stderr,"rare m=%d done nodes=%lld t=%.0f\n",m,nodes,now()-t0);
    }
    for(int x=0;x<100;x++){req[x]=-1;forb[x]=0;}
  }
  fprintf(out,"DONE nodes=%lld leaves=%lld fams=%lld survivors=%.0f frac=%.9g\n",nodes,leaves,fams,famcount,frac);
  return 0;
}
