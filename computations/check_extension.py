"""Independent squared-argument arithmetic and primitive; 90 digit intervals."""
import importlib.util,json,gzip,time,os,sys,resource
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode=True
p=Path('computations/restriction_v4/audit.py')
spec=importlib.util.spec_from_file_location('other',p);a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)
start=time.monotonic();p=Path('computations/restriction_extension')
row=json.loads((p/'proposal.json').read_text());cert=json.loads((p/'block_certificate.json').read_text());n=row['L'];U=row['U']
assert (n,U)==(cert['L'],cert['U'])
beta,nw,slope,intercept,cost=a.params(row,n);alpha=a.I(cert['alpha'])
assert a.hi(a.pair(F(1,n))[1]+slope)<a.lo(alpha+sum((w for t,w in nw),a.I(0)))
with gzip.open(p/cert['covering'],'rt') as f:cover=json.load(f)
cursor=F(1,U);minimum=None
for left,right in cover:
 l,r=F(left),F(right);assert l==cursor and l<r<=1
 lhs=a.pair(l)[0]+a.shift(l,r,n)
 rhs=alpha-beta*a.I(r)+sum((w*a.iv.exp(-t*a.I(r)) for t,w in nw),a.I(0))
 slack=a.lo(rhs)-a.hi(lhs);assert slack>0,(l,r);minimum=min(slack,minimum) if minimum is not None else slack;cursor=r
assert cursor==1 and len(cover)==cert['cells']
bound=(alpha+intercept-beta*a.I(F(n-1,2*n))+cost)/a.iv.ln(2)
assert a.hi(bound)<F(cert['rational_upper'])
# Independent termwise interval evaluation of Proposition 4.5.
tail=json.loads((p/'tail_certificate.json').read_text());n=tail['n0'];av=F(tail['a']);assert 0<av<=F(1,5) and av*n>=40
I,iv=a.I,a.iv;N=I(n-1);v=I(av);r=I(F(23,48));c=N/n*(1-iv.exp(-2*v))/(2*v)
vT=N*(1-iv.exp(-2*v))/4;v20=N*(1-iv.exp(-I(F(40,n))))/4
b=iv.sqrt(2/(iv.pi*N));z=40*r*r;A=-I(F(1,2))-iv.ln(v)+v
terms=[b*(3+iv.ln(1+1/(b*b))+2*iv.ln(2)),(6+4*iv.ln(n))/n,I(F(14,5)),I(F(273,100*n)),I(F(48,23))/iv.sqrt(40*iv.pi)+iv.exp(-z)/z,a.phi(20*c)*iv.ln(1/c),iv.exp(v)*A*(I(F(24,23))/iv.sqrt(2*iv.pi*vT)+iv.exp(-2*r*r*vT)),I(F(89,50))*a.phi(v20)]
for t,old in zip(terms,tail['terms_nats']):assert F(old[0])<=a.lo(t)<=a.hi(t)<=F(old[1])
tb=sum(terms,I(0))/iv.ln(2);assert a.hi(tb)<F(tail['rational_upper'])
# Integral quadrature checks the special-function values by another method.
import mpmath as mp
mp.mp.dps=85
M=lambda x:mp.mpf(F(x).numerator)/F(x).denominator
n=tail['n0'];N=mp.mpf(n-1);v=M(av);r=M(F(23,48));c=N/n*(-mp.expm1(-2*v))/(2*v)
vT=N*(-mp.expm1(-2*v))/4;v20=N*(-mp.expm1(-mp.mpf(40)/n))/4
phi=lambda x:2/mp.pi*mp.quad(lambda t:mp.exp(-2*x*mp.sin(t)**2),[0,mp.pi/2])
b=mp.sqrt(2/(mp.pi*N));z=40*r*r;A=-mp.mpf(1)/2-mp.log(v)+v
qt=(b*(3+mp.log1p(1/b**2)+2*mp.log(2))+(6+4*mp.log(n))/n+M(F(14,5))+M(F(273,100*n))+48/(23*mp.sqrt(40*mp.pi))+mp.exp(-z)/z+phi(20*c)*mp.log(1/c)+mp.exp(v)*A*(24/(23*mp.sqrt(2*mp.pi*vT))+mp.exp(-2*r*r*vT))+M(F(89,50))*phi(v20))/mp.log(2)
assert M(tail['bound_interval'][0])<=qt<=M(tail['bound_interval'][1])
out={'passed':True,'cells':len(cover),'interval_digits':iv.dps,'minimum_cell_slack_nats':str(minimum),'block_bound_upper':str(a.hi(bound)),'tail_bound_upper':str(a.hi(tb)),'tail_quadrature':str(qt),'seconds':time.monotonic()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'nice':os.getpriority(os.PRIO_PROCESS,0)}
(p/'audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
