"""Independent big-integer audit, exact prime-log signatures, and Walsh checks."""
import sys, importlib.util, json, math
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
sys.dont_write_bytecode=True
src=Path(__file__).resolve().parent/'reference/small_dimensions.py'
spec=importlib.util.spec_from_file_location('ref32',src); ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
@lru_cache(None)
def matrix(n): return ref.kraw(n)
def spectrum(t):return ref.transform([1-2*int(x) for x in t],matrix(len(t)-1))
def data(t):
 n=len(t)-1;A=spectrum(t);c=[math.comb(n,s) for s in range(n+1)]
 e=ref.ent_bounds(n,A,c)
 return n,A,c,e

def block(t,k):
 n,A,c,e=data(t);lo,hi=e
 coeff={2:2*k*4**n}
 def add(A,c,scale):
  for a,m in zip(A,c):
   for p,v in ref.factors(abs(a)) if a else ():
    coeff[p]=coeff.get(p,0)+scale*m*a*a*v
 add(A,c,-2)
 for r in range(k+1):
  _,B,d,(l,h)=data(t[r:r+n-k+1]);z=2**k*math.comb(k,r)
  lo-=z*h;hi-=z*l;add(B,d,2*z)
 edges=sum(math.comb(n-1,w) for w in range(n) if t[w]!=t[w+1])
 den=k*edges*2**(n+1)
 sig={str(p):str(Fraction(v,den)) for p,v in sorted(coeff.items()) if v}
 return Fraction(lo,ref.Q*den),Fraction(hi,ref.Q*den),sig
rows=[json.loads(x) for x in open('exhaustive.jsonl')]
for row in rows:
 n,k=row['n'],row['k'];sig0=None
 for t in row['maximizers_mod_output_complement']:
  lo,hi,sig=block(t,k)
  if sig0 is None:
   sig0=sig
   row['certified_ratio_lower']=str(lo);row['certified_ratio_upper']=str(hi)
   row['decimal_enclosure']=[ref.decimal_bound(lo.numerator,lo.denominator,18),ref.decimal_bound(hi.numerator,hi.denominator,18,True)]
   assert Fraction(int(row['ratio_lo_numerator']),int(row['denominator']))<=lo
   assert Fraction(int(row['ratio_hi_numerator']),int(row['denominator']))>=hi
  assert sig==sig0
 row['exact_ratio_prime_log_coefficients']=sig0
 row['all_maximizers'] = sorted(set(row['maximizers_mod_output_complement']+[''.join('1' if x=='0' else '0' for x in t) for t in row['maximizers_mod_output_complement']]))
 row['strict_exclusion_and_symbolic_ties_verified']=True
Path('certified_maxima.json').write_text(json.dumps(rows,indent=2)+'\n')
# Entire Boolean cubes, independent of Krawtchouk, for small dimensions.
def walsh(t):
 n=len(t)-1;a=[1-2*int(t[x.bit_count()]) for x in range(1<<n)];d=1
 while d<len(a):
  for s in range(0,len(a),2*d):
   for j in range(d):u,v=a[s+j],a[s+j+d];a[s+j],a[s+j+d]=u+v,u-v
  d*=2
 assert sum(x*x for x in a)==4**n
 for x,v in enumerate(a):assert v==spectrum(t)[x.bit_count()]
 return sum((x*x/4**n)*math.log2(4**n/(x*x)) for x in a if x)
checks=0
for n in range(1,7):
 maxima={k:-1 for k in range(2,min(5,n)+1)}
 for mask in range(1,2**n):
  t=''.join(str((mask>>w)&1) for w in range(n+1));h=walsh(t)
  edges=sum(math.comb(n-1,w) for w in range(n) if t[w]!=t[w+1]);i1=edges/2**(n-1)
  for k in maxima:
   delta=h-sum(math.comb(k,r)*walsh(t[r:r+n-k+1]) for r in range(k+1))/2**k
   maxima[k]=max(maxima[k],delta/(k*i1))
  checks+=1
 for k,v in maxima.items():
  row=next(r for r in rows if r['n']==n and r['k']==k)
  assert abs(v-float(Fraction(row['certified_ratio_lower'])))<1e-12
# OR closed form independently checks every global finite-range winner.
def orh(n):
 if n==0:return 0.
 q=2.**(1-n);p0=(1-q)**2
 return -p0*math.log2(p0) if n==1 else (-p0*math.log2(p0)+(2**n-1)*q*q*(2*n-2))
for k in range(2,6):
 n=k+2; row=next(r for r in rows if r['n']==n and r['k']==k)
 q=(orh(n)-orh(n-k)/2**k)/(k*2.**(1-n))
 assert abs(q-float(Fraction(row['certified_ratio_lower'])))<1e-12
Path('audit.json').write_text(json.dumps({'rows_audited':len(rows),'independent_Walsh_nonconstant_tables':checks,'maximum_dimension':20,'interval_bits_cpp':56,'interval_bits_python':ref.BITS,'atanh_terms_cpp':30,'atanh_terms_python':ref.TERMS,'all_checks_passed':True},indent=2)+'\n')
with open('maxima.csv','w') as f:
 f.write('k,n,ratio_lower,ratio_upper,maximizers_mod_output_complement\n')
 for r in rows:f.write(f"{r['k']},{r['n']},{r['decimal_enclosure'][0]},{r['decimal_enclosure'][1]},{'|'.join(r['maximizers_mod_output_complement'])}\n")
print(Path('audit.json').read_text())
