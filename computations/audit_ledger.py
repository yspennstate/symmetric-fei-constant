"""Independent endpoint and Decimal reconstruction of the new ledger."""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal,localcontext
import csv,gzip,json,hashlib,time,os
start=time.monotonic();F=Fraction
b=Path('computations/restriction_v4');out=json.loads(Path('computations/restriction_extension/ledger_certificate.json').read_text())
old=json.loads((b/'ledger_certificate.json').read_text())
raw=json.loads(Path('computations/block_census_24/certified_maxima.json').read_text())
assert len(raw)==86
# Reconstruct the seed as a shortest path in integer-valued n*C_n, avoiding
# the producer's recurrence directly in C_n.
scaled=[F(0),F(0)]
for n in range(2,25):
 candidates=[scaled[n-r['k']]+r['k']*F(r['certified_ratio_upper']) for r in raw if r['n']==n]
 scaled.append(min(candidates))
 assert scaled[n]/n==F(out['small_dimensions'][str(n)])
assert scaled[24]==24*F(out['base_ratio'])
# Get bounds afresh from original certificate files rather than copying the
# producer's assembled segment list.
finite=json.loads((b/'finite_certificate.json').read_text())['dimensions']
blocks=json.loads((b/'blocks_certificate.json').read_text())['blocks']
ext=json.loads(Path('computations/restriction_extension/block_certificate.json').read_text())
segments=[(r['n'],r['n'],F(r['rational_upper'])) for r in finite if r['n']>24]
segments += [(r['L'],r['U'],F(r['rational_upper'])) for r in blocks]+[(ext['L'],ext['U'],F(ext['rational_upper']))]
acc=scaled[24];peak=max((scaled[n]/n,n) for n in range(1,25));cursor=25
for i,(L,U,B) in enumerate(segments):
 assert L==cursor;cursor=U+1
 # On a constant segment a(n)=B+(acc-(L-1)B)/n.
 endpoints=[(acc/(L-1),L-1),((acc+(U-L+1)*B)/U,U)]
 peak=max(peak,*endpoints);acc+=(U-L+1)*B
 assert str(acc)==out['segments'][i]['cumulative_sum_end']
assert str(peak[0])==out['finite_maximum_rational'] and peak[1]==out['binding_dimension']
tail=json.loads(Path('computations/restriction_extension/tail_certificate.json').read_text());assert cursor-1==tail['n0']==out['cutoff']
maximum=max(peak[0],F(tail['rational_upper']));assert str(maximum)==out['all_dimension_constant_rational']
safe=F(out['safe_rational']);assert maximum<=safe and safe-F(1,10**6)<maximum
# Validate every emitted head row independently against exact block sums.
head_checks=0
with gzip.open('computations/restriction_extension/head_ledger.csv.gz','rt') as f:
 for row in csv.DictReader(f):
  n=int(row['dimension']);S=scaled[24]
  for L,U,B in segments:
   if L>n:break
   S+=max(0,min(U,n)-L+1)*B
  assert S==F(row['cumulative_upper']) and S/n==F(row['H_over_I_upper']);head_checks+=1
with localcontext() as ctx:
 ctx.prec=100
 def D(x):return Decimal(x.numerator)/Decimal(x.denominator)
 S=D(scaled[24]);n=peak[1]
 for L,U,B in segments:
  if L>n:break
  S+=max(0,min(U,n)-L+1)*D(B)
 decimal_maximum=S/n
 assert abs(decimal_maximum-D(peak[0]))<Decimal('1e-90')
 cutoff_decimal=D(scaled[24])
 for L,U,B in segments:cutoff_decimal+=(U-L+1)*D(B)
 cutoff_decimal/=out['cutoff']
 assert abs(cutoff_decimal-D(acc/out['cutoff']))<Decimal('1e-90')
# Every row in independent new-census replay matches the exclusion counts.
new={(r['n'],r['k']):r for r in raw if r['n']>=21}
rechecks=[json.loads(s) for s in Path('computations/block_census_24/exclusion_recheck.jsonl').read_text().splitlines()]
assert len(rechecks)==16
for r in rechecks:
 orig=new[r['n'],r['k']]
 assert r['strict_exclusions']==orig['strictly_excluded_mod_complement'] and r['ties']==orig['overlapping_mod_complement']
result={'passed':True,'independent_method':'integer n*C_n shortest paths, rational segment endpoints, 100-digit decimal binding value, full head-row verification','cutoff_decimal':str(cutoff_decimal),'head_rows_checked':head_checks,'new_exclusion_certificates_checked':len(rechecks),'decimal_maximum':str(decimal_maximum),'finite_maximum_rational':str(peak[0]),'binding_dimension':peak[1],'tail_gap_below_finite_maximum':str(peak[0]-F(tail['rational_upper'])),'seconds':time.monotonic()-start,'nice':os.getpriority(os.PRIO_PROCESS,0)}
Path('computations/restriction_extension/ledger_audit.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
