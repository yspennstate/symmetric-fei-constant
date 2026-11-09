"""Exact Walsh certificates for the two general-function restriction examples."""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import json
from mpmath import iv
iv.dps=80

def walsh(a):
 a=a.copy();d=1
 while d<len(a):
  for start in range(0,len(a),2*d):
   for j in range(start,start+d):
    u,v=a[j],a[j+d];a[j],a[j+d]=u+v,u-v
  d*=2
 return a

def factors(a):
 out={};p=2
 while p*p<=a:
  while a%p==0:out[p]=out.get(p,0)+1;a//=p
  p+=1
 if a>1:out[a]=out.get(a,0)+1
 return out

def entropy(a):
 n=(len(a)).bit_length()-1;A=walsh(a)
 assert sum(z*z for z in A)==4**n
 sig=Counter({2:F(2*n)})
 for z,count in Counter(abs(z) for z in A).items():
  if z:
   for p,k in factors(z).items():sig[p]-=F(2*count*z*z*k,4**n)
 return dict(sig),A

def value(sig):
 out=iv.mpf(0)
 for p,q in sig.items():out+=(iv.mpf(q.numerator)/q.denominator)*iv.ln(p)/iv.ln(2)
 return out

def restrict(a,i,b):
 return [a[(x&((1<<i)-1))|((x>>i)<<(i+1))|(b<<i)] for x in range(len(a)//2)]

def endpoints(x):
 out=[]
 for sign,m,e,_ in x._mpi_:
  q=F(int(m))*F(2)**int(e);out.append(str(-q if sign else q))
 return out

rows=[]
for n,name,threshold in [(6,'majority of three AND pairs',2),(16,'OR of eight AND pairs',1)]:
 a=[1-2*int(sum(((x>>(2*j))&3)==3 for j in range(n//2))>=threshold) for x in range(1<<n)]
 sig,A=entropy(a);common=None;influences=[]
 for i in range(n):
  count=sum(a[x]!=a[x^(1<<i)] for x in range(1<<n))
  influence=F(count,1<<n)
  assert influence==F(sum(z*z for j,z in enumerate(A) if j&(1<<i)),4**n)
  delta=Counter(sig)
  for b in [0,1]:
   sub,_=entropy(restrict(a,i,b))
   for p,q in sub.items():delta[p]-=q/2
  ratio={p:q/influence for p,q in delta.items() if q}
  if common is None:common=ratio
  assert ratio==common
  influences.append(str(influence))
 v=value(common)
 expected=F('5.1948' if n==6 else '6.0037')
 assert v.a>iv.mpf(expected.numerator)/expected.denominator
 rows.append({'variables':n,'function':name,'all_coordinate_ratios_equal_exactly':True,
 'coordinate_influences':influences,'ratio_prime_log_coefficients':{str(p):str(q) for p,q in sorted(common.items())},
 'ratio_interval_rational':endpoints(v),'ratio_interval_decimal':str(v),'interval_digits':iv.dps})
 print(name,str(v))
Path(__file__).with_name('examples_v4.json').write_text(json.dumps(rows,indent=2)+'\n')
