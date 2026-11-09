// Exact integer Krawtchouk enumeration and fixed-point rational logarithm enclosures.
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
using U=__uint128_t; using Z=__int128_t;
const uint64_t Q=1ULL<<56;
struct Iv { Z lo=0,hi=0; };
struct Log {uint64_t lo=0,hi=0;};
U ceildiv(U a,U b){return (a+b-1)/b;}
uint64_t mullo(uint64_t a,uint64_t b){return (U(a)*b)>>56;}
uint64_t mulhi(uint64_t a,uint64_t b){return ceildiv(U(a)*b,Q);}
Log lnrat(uint64_t a,uint64_t b){
 // z=a/b<=1/3; ln((1+z)/(1-z))=2 sum z^(2j+1)/(2j+1).
 assert(3*a<=b);
 if(!a)return {};
 uint64_t zl=U(Q)*a/b,zu=ceildiv(U(Q)*a,b);
 uint64_t l2=mullo(zl,zl),u2=mulhi(zu,zu),pl=zl,pu=zu,sl=0,su=0;
 for(int j=0;j<30;j++){
   sl+=pl/(2*j+1); su+=uint64_t(ceildiv(pu,2*j+1));
   pl=mullo(pl,l2);pu=mulhi(pu,u2);
 }
 // Tail <=9/(4*61*3^61)<1/Q. All intermediates outward rounded.
 return {2*sl,2*su+1};
}
vector<Log> logs;
Log lg(int a){return logs[a];}
string str(Z a){if(!a)return "0";bool neg=a<0;if(neg)a=-a;string s;while(a){s+=char('0'+a%10);a/=10;}if(neg)s+='-';reverse(s.begin(),s.end());return s;}
string tab(uint32_t mask,int n){string s;for(int w=0;w<=n;w++)s+=char('0'+((mask>>w)&1));return s;}
uint64_t C[21][21];
vector<vector<int64_t>> matrix(int n){
 vector<vector<int64_t>> a(n+1,vector<int64_t>(n+1));
 for(int s=0;s<=n;s++)for(int w=0;w<=n;w++)for(int j=max(0,w-(n-s));j<=min(s,w);j++)
 a[s][w]+=(j%2?-1:1)*int64_t(C[s][j]*C[n-s][w-j]);
 return a;
}
int main(int argc,char**argv){
 int maxn=argc>1?stoi(argv[1]):20;
 assert(maxn<=20);
 U three61=1;for(int j=0;j<61;j++)three61*=3;assert(U(9)*Q<U(4)*61*three61);
 for(int n=0;n<=20;n++){C[n][0]=C[n][n]=1;for(int j=1;j<n;j++)C[n][j]=C[n-1][j-1]+C[n-1][j];}
 auto ln2=lnrat(1,3);logs.resize((1<<maxn)+1);
 for(int a=1;a<=1<<maxn;a++){
  int b=31-__builtin_clz(a),p=1<<b;auto z=lnrat(a-p,a+p);
  logs[a]={uint64_t(b)*Q+uint64_t(U(Q)*z.lo/ln2.hi),uint64_t(b)*Q+uint64_t(ceildiv(U(Q)*z.hi,ln2.lo))};
 }
 vector<vector<Iv>> ent(maxn+1);ent[0].resize(1);
 ofstream out("exhaustive.jsonl");out<<setprecision(20);
 for(int n=1;n<=maxn;n++){
  auto start=chrono::steady_clock::now();auto M=matrix(n);int count=1<<n;ent[n].resize(count);
  vector<int64_t>A(n+1);A[0]=1<<n;
  uint32_t prev=0;vector<int>F(n+1,1);
  vector<Iv>best(6);vector<uint64_t>be(6);vector<vector<uint32_t>>ties(6),overlaps(6);
  for(int ctr=0;ctr<count;ctr++){
   uint32_t mask=ctr^(ctr>>1);
   if(ctr){int w=__builtin_ctz(mask^prev);int d=-2*F[w];F[w]=-F[w];for(int s=0;s<=n;s++)A[s]+=d*M[s][w];}prev=mask;
   Iv e;U parse=0;uint64_t jfour=0,edges=0;
   for(int w=0;w<n;w++)if(F[w]!=F[w+1])edges+=C[n-1][w];
   for(int s=0;s<=n;s++){
    uint64_t a=llabs(A[s]),mass=C[n][s]*a*a;parse+=mass;
    if(s)jfour+=C[n-1][s-1]*a*a;
    if(a){auto z=lg(a);e.lo+=Z(mass)*(Z(2*n)*Q-2*z.hi);e.hi+=Z(mass)*(Z(2*n)*Q-2*z.lo);}
   }
   assert(parse==(U(1)<<(2*n)));assert(jfour==(edges<<(n+1)));ent[n][mask]=e;
   if(!edges)continue;
   for(int k=2;k<=min(5,n);k++){
    Iv delta=e;int m=n-k;uint32_t all=(1U<<(m+1))-1;
    for(int r=0;r<=k;r++){
     uint32_t sub=(mask>>r)&all;sub=min(sub,sub^all);Iv z=ent[m][sub];
     delta.lo-=Z(1<<k)*C[k][r]*z.hi;delta.hi-=Z(1<<k)*C[k][r]*z.lo;
    }
    if(!be[k]||delta.lo*be[k]>best[k].hi*edges){best[k]=delta;be[k]=edges;ties[k]={mask};overlaps[k].clear();}
    else if(delta.hi*be[k]>=best[k].lo*edges){ties[k].push_back(mask);}
   }
  }
  // A second complete pass compares every table to the final reference.
  vector<uint64_t>excluded(6),overlap(6);
  for(uint32_t mask=1;mask<uint32_t(count);mask++){
   uint64_t edges=0;for(int w=0;w<n;w++)if(((mask>>w)^(mask>>(w+1)))&1)edges+=C[n-1][w];
   assert(edges);
   for(int k=2;k<=min(5,n);k++){
    Iv delta=ent[n][mask];int m=n-k;uint32_t all=(1U<<(m+1))-1;
    for(int r=0;r<=k;r++){
     uint32_t sub=(mask>>r)&all;sub=min(sub,sub^all);auto z=ent[m][sub];
     delta.lo-=Z(1<<k)*C[k][r]*z.hi;delta.hi-=Z(1<<k)*C[k][r]*z.lo;
    }
    if(delta.hi*be[k]<best[k].lo*edges)excluded[k]++;
    else {assert(find(ties[k].begin(),ties[k].end(),mask)!=ties[k].end());overlap[k]++;}
   }
  }
  for(int k=2;k<=min(5,n);k++){
   assert(excluded[k]+overlap[k]==uint64_t(count-1));
   // All possible overlapping maxima are exported for exact symbolic tie checking.
   Z den=Z(Q)*k*be[k]*(1ULL<<(n+1));
   out<<"{\"strictly_excluded_mod_complement\":"<<excluded[k]<<",\"overlapping_mod_complement\":"<<overlap[k]<<",\"n\":"<<n<<",\"k\":"<<k<<",\"tables_covered\":"<<(1ULL<<(n+1))<<",\"enumerated_up_to_output_complement\":"<<count<<",\"ratio_lo_numerator\":\""<<str(best[k].lo)<<"\",\"ratio_hi_numerator\":\""<<str(best[k].hi)<<"\",\"denominator\":\""<<str(den)<<"\",\"ratio_approx\":"<<(long double)best[k].lo/(long double)den<<",\"maximizers_mod_output_complement\":[";
   for(size_t i=0;i<ties[k].size();i++){if(i)out<<",";out<<"\""<<tab(ties[k][i],n)<<"\"";}out<<"],\"seconds_for_dimension\":"<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<"}\n";out.flush();
  }
  cerr<<"completed n="<<n<<" canonical tables="<<count<<" seconds="<<chrono::duration<double>(chrono::steady_clock::now()-start).count()<<endl;
 }
}
