// Independent census replay: ordinary binary traversal, 48-bit logarithms,
// checked 128-bit comparisons, and a prescribed reference for each new pair.
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <vector>
#include <sys/resource.h>
using namespace std;
using u128=__uint128_t;using i128=__int128_t;
const uint64_t SCALE=1ULL<<48;
struct Bounds{i128 lower=0,upper=0;};
struct Logs{uint64_t lower,upper;};
#include "references.h"
uint64_t binomial[25][25];
uint64_t up(u128 a,uint64_t b){return a/b+(a%b!=0);}
Logs ln_unit(uint64_t a,uint64_t b){
 assert(3*a<=b);if(a==0)return {0,0};
 // 2 atanh(a/b), with outward intervals for every monomial. 40 terms,
 // remainder < 9/(4*81*3^81) < 1/SCALE.
 uint64_t l=u128(a)*SCALE/b,h=up(u128(a)*SCALE,b);
 uint64_t l2=u128(l)*l/SCALE,h2=up(u128(h)*h,SCALE),pl=l,ph=h,sl=0,sh=0;
 for(int j=0;j<40;j++){
  sl+=pl/(2*j+1);sh+=up(ph,2*j+1);
  pl=u128(pl)*l2/SCALE;ph=up(u128(ph)*h2,SCALE);
 }
 return {2*sl,2*sh+1};
}
vector<Logs> logs;
vector<vector<int64_t>> krawtchouk(int n){
 vector<vector<int64_t>> m(n+1,vector<int64_t>(n+1));
 // Generating polynomials (1-z)^s (1+z)^(n-s), multiplied directly.
 for(int s=0;s<=n;s++){
  vector<int64_t> p(n+1);p[0]=1;
  for(int j=0;j<n;j++)for(int w=j+1;w>=1;w--)p[w]+=(j<s?-1:1)*p[w-1];
  m[s]=p;
 }
 return m;
}
Bounds entropy(int n,const vector<int64_t>& a){
 Bounds z;u128 parseval=0;
 for(int s=0;s<=n;s++){
  uint64_t x=a[s]<0?-a[s]:a[s];u128 mass=u128(binomial[n][s])*x*x;parseval+=mass;
  if(x){z.lower+=i128(mass)*(i128(2*n)*SCALE-2*logs[x].upper);z.upper+=i128(mass)*(i128(2*n)*SCALE-2*logs[x].lower);}
 }
 assert(parseval==(u128(1)<<(2*n)));return z;
}
i128 multiply(i128 a,uint64_t b){i128 z;assert(!__builtin_mul_overflow(a,i128(b),&z));return z;}
int main(){
 const int N=24;auto start=chrono::steady_clock::now();
 for(int n=0;n<=N;n++)for(int j=0;j<=n;j++)binomial[n][j]=(j==0||j==n)?1:binomial[n-1][j-1]+binomial[n-1][j];
 // The tail inequality is checked with integer arithmetic via 3^60 already:
 // 9*SCALE < 4*81*3^81 follows from the stronger 9*SCALE < 4*81*3^60.
 u128 th=1;for(int j=0;j<60;j++)th*=3;assert(u128(9)*SCALE<u128(4)*81*th);
 Logs l2=ln_unit(1,3);logs.resize((1<<N)+1);
 for(uint32_t x=1;x<=uint32_t(1<<N);x++){
  int b=31-__builtin_clz(x);uint64_t p=1ULL<<b;Logs v=ln_unit(x-p,x+p);
  logs[x]={b*SCALE+uint64_t(u128(v.lower)*SCALE/l2.upper),b*SCALE+up(u128(v.upper)*SCALE,l2.lower)};
 }
 vector<vector<Bounds>> ent(N+1);ent[0].resize(1);
 ofstream out("exclusion_recheck.jsonl");uint64_t total=0;
 for(int n=1;n<=N;n++){
  auto mat=krawtchouk(n);uint32_t count=1U<<n;ent[n].resize(count);
  array<Bounds,6> ref;array<uint64_t,6> refedges={};array<uint64_t,6> excluded={},ties={};
  if(n>=21)for(int k=2;k<=5;k++){
   uint32_t mask=REFERENCE[n][k];vector<int64_t>a(n+1);
   for(int s=0;s<=n;s++)for(int w=0;w<=n;w++)a[s]+=(((mask>>w)&1)?-1:1)*mat[s][w];
   ref[k]=entropy(n,a);
   for(int w=0;w<n;w++)if(((mask>>w)^(mask>>(w+1)))&1)refedges[k]+=binomial[n-1][w];
   for(int r=0;r<=k;r++){
    int m=n-k;uint32_t all=(1U<<(m+1))-1,sub=(mask>>r)&all;sub=min(sub,sub^all);
    auto e=ent[m][sub];ref[k].lower-=i128(1U<<k)*binomial[k][r]*e.upper;ref[k].upper-=i128(1U<<k)*binomial[k][r]*e.lower;
   }
  }
  vector<int64_t>a(n+1);a[0]=1U<<n;uint32_t prev=0;
  for(uint32_t mask=0;mask<count;mask++){
   // Binary carries can flip multiple coordinates. This visits the tables in
   // a different order than the producer's one-bit Gray transitions.
   uint32_t change=mask^prev;
   while(change){int w=__builtin_ctz(change);change&=change-1;int delta=((mask>>w)&1)?-2:2;for(int s=0;s<=n;s++)a[s]+=delta*mat[s][w];}
   prev=mask;auto e=entropy(n,a);ent[n][mask]=e;
   if(n<21||mask==0)continue;
   uint64_t edges=0;for(int w=0;w<n;w++)if(((mask>>w)^(mask>>(w+1)))&1)edges+=binomial[n-1][w];assert(edges);
   for(int k=2;k<=5;k++){
    Bounds delta=e;
    for(int r=0;r<=k;r++){
     int m=n-k;uint32_t all=(1U<<(m+1))-1,sub=(mask>>r)&all;sub=min(sub,sub^all);auto z=ent[m][sub];
     delta.lower-=i128(1U<<k)*binomial[k][r]*z.upper;delta.upper-=i128(1U<<k)*binomial[k][r]*z.lower;
    }
    if(multiply(delta.upper,refedges[k])<multiply(ref[k].lower,edges))excluded[k]++;
    else{
     // Every nonexcluded table is an audited exact symbolic tie.
     uint32_t refl=REFERENCE[n][k],reversed=0;
     for(int w=0;w<=n;w++)reversed|=((refl>>w)&1)<<(n-w);
     uint32_t all=(1U<<(n+1))-1;reversed=min(reversed,reversed^all);
     assert(mask==refl||mask==reversed);ties[k]++;
    }
   }
   total++;
  }
  if(n>=21)for(int k=2;k<=5;k++){
   assert(excluded[k]+ties[k]==count-1);assert(ties[k]==2);
   out<<"{\"n\":"<<n<<",\"k\":"<<k<<",\"strict_exclusions\":"<<excluded[k]<<",\"ties\":"<<ties[k]<<"}\n";out.flush();
  }
  cerr<<"rechecked n="<<n<<endl;
 }
 struct rusage usage;getrusage(RUSAGE_SELF,&usage);
 ofstream metrics("exclusion_metrics.json");metrics<<"{\"new_nonconstant_tables\":"<<total<<",\"fixed_point_bits\":48,\"seconds\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<",\"peak_rss_bytes\":"<<usage.ru_maxrss<<",\"nice\":"<<getpriority(PRIO_PROCESS,0)<<"}\n";
}
