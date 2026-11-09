"""Full per-dimension Fraction ledger, with two exact segment cross-checks."""
import os,json,gzip,csv,time,resource,hashlib
from fractions import Fraction as F
from decimal import Decimal,localcontext
from pathlib import Path
start=time.monotonic();base=Path('computations/restriction_v4')
old=json.loads((base/'ledger_certificate.json').read_text())
rows=json.loads(Path('computations/block_census_24/certified_maxima.json').read_text());B={(r['k'],r['n']):F(r['certified_ratio_upper']) for r in rows}
C={0:F(0),1:F(0)};choices={}
for n in range(2,25):C[n],choices[n]=min(((k*B[k,n]+(n-k)*C[n-k])/n,k) for k in range(2,min(n,5)+1))
for n in range(21):assert C[n]==F(old['small_dimensions'][str(n)])
segments=[{'L':max(25,s['L']),'U':s['U'],'upper':s['upper'],'method':s['method']} for s in old['segments'] if s['U']>=25]
ex=json.loads(Path('computations/restriction_extension/block_certificate.json').read_text());segments.append({'L':ex['L'],'U':ex['U'],'upper':ex['rational_upper'],'method':'new continuous extension'})
tail=json.loads(Path('computations/restriction_extension/tail_certificate.json').read_text());assert segments[-1]['U']==tail['n0']
total=24*C[24];best=(max(C.values()),max(C,key=C.get));cursor=25
with gzip.open('computations/restriction_extension/head_ledger.csv.gz','wt') as f:
 w=csv.writer(f);w.writerow(['dimension','one_step_upper','cumulative_upper','H_over_I_upper'])
 for s in segments:
  assert s['L']==cursor;cursor=s['U']+1;u=F(s['upper']);before=total
  endpoint_best=max((before/(s['L']-1),s['L']-1),((before+(s['U']-s['L']+1)*u)/s['U'],s['U']))
  local_best=(before/(s['L']-1),s['L']-1)
  for n in range(s['L'],s['U']+1):
   total+=u;q=total/n;best=max(best,(q,n));local_best=max(local_best,(q,n))
   if n<=16384:w.writerow([n,str(u),str(total),str(q)])
  assert total==before+(s['U']-s['L']+1)*u and local_best==endpoint_best
  s['cumulative_sum_end']=str(total)
  if s['U']>=1048576:print('ledger reached',s['U'],flush=True)
constant=max(best[0],F(tail['rational_upper']));ceiling=(constant.numerator*10**6+constant.denominator-1)//constant.denominator
safe=F(ceiling,10**6)
def dec(q):
 with localcontext() as c:c.prec=85;return str(Decimal(q.numerator)/Decimal(q.denominator))
result={'all_dimension_constant_rational':str(constant),'all_dimension_constant_decimal':dec(constant),'safe_rational':str(safe),'safe_decimal':dec(safe),'binding_dimension':best[1],'finite_maximum_rational':str(best[0]),'finite_maximum_decimal':dec(best[0]),'base_dimension':24,'base_ratio':str(C[24]),'small_dimensions':{str(n):str(v) for n,v in C.items()},'small_block_choices':choices,'tail_upper':tail['rational_upper'],'cutoff':tail['n0'],'ledger_steps':tail['n0']-24,'segments':segments,'seconds':time.monotonic()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'nice':os.getpriority(os.PRIO_PROCESS,0)}
result['improvement_over_published']=dec(F('4.657435')-safe);result['improvement_over_old_exact_maximum']=dec(F(old['all_dimension_constant_rational'])-constant)
inputs=['computations/block_census_24/certified_maxima.json','computations/restriction_extension/block_certificate.json','computations/restriction_extension/tail_certificate.json','computations/restriction_v4/ledger_certificate.json']
result['input_sha256']={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in inputs}
Path('computations/restriction_extension/ledger_certificate.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ('segments','small_dimensions')},indent=2),flush=True)
