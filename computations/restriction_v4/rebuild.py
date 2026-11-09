"""Rebuild the complete exact restriction ledger from local certificates.

Run with --verify to replay all interval certificates first. The ledger loop
adds every per-d bound as a Fraction, even in constant-valued block ranges.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import argparse,csv,gzip,hashlib,json,resource,subprocess,sys,time
from decimal import Decimal,localcontext
from fractions import Fraction as F
from pathlib import Path

def dec(q):
    with localcontext() as c:
        c.prec=75;return str(Decimal(q.numerator)/Decimal(q.denominator))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    start=time.monotonic()
    if args.verify:subprocess.run([sys.executable,'certify.py','--replay'],check=True)
    census=json.loads(Path('source_inputs/small_maxima.json').read_text())
    K={(r['k'],r['n']):F(r['certified_ratio_upper']) for r in census}
    C={0:F(0),1:F(0)};small_steps={}
    for n in range(2,21):
        C[n],small_steps[n]=min(((k*K[k,n]+(n-k)*C[n-k])/n,k) for k in range(2,min(5,n)+1))
    assert C[20]<=F(json.loads(Path('source_inputs/source_chain.json').read_text())['base_H_over_I_upper'])
    exact=json.loads(Path('finite_certificate.json').read_text())
    assert exact['coefficients_checked']==32520
    dims=exact['dimensions'];assert [r['n'] for r in dims]==list(range(16,256))
    blocks=json.loads(Path('blocks_certificate.json').read_text());assert blocks['complete']
    segments=[{'L':r['n'],'U':r['n'],'upper':r['rational_upper'],'method':'exact dimension'} for r in dims if r['n']>=21]
    segments += [{'L':r['L'],'U':r['U'],'upper':r['rational_upper'],'method':'continuous block'} for r in blocks['blocks']]
    cursor=21
    for r in segments:assert r['L']==cursor and r['U']>=r['L'];cursor=r['U']+1
    cutoff=cursor-1;tail=json.loads(Path('tail_certificate.json').read_text());assert tail['n0']==cutoff
    total=20*C[20];best=(total/20,20);checkpoints=[]
    with gzip.open('head_ledger.csv.gz','wt') as f:
        writer=csv.writer(f);writer.writerow(['dimension','U_d_upper','cumulative_nats_free_bits_bound','H_over_I_upper'])
        for seg in segments:
            bound=F(seg['upper']);before=total
            blockbest=(total/(seg['L']-1),seg['L']-1)
            for n in range(seg['L'],seg['U']+1):
                total+=bound;q=total/n;best=max(best,(q,n));blockbest=max(blockbest,(q,n))
                if n<=16384:writer.writerow([n,str(bound),str(total),str(q)])
            # Independent exact multiplication checks all additions in this segment.
            assert total==before+(seg['U']-seg['L']+1)*bound
            # A constant-bound segment's average is monotone, so its maximum is
            # attained at an endpoint. Check the enumerated maximum a second way.
            assert blockbest==max((before/(seg['L']-1),seg['L']-1),(total/seg['U'],seg['U']))
            seg['cumulative_sum_end']=str(total)
            if seg['method']=='continuous block':
                checkpoints.append({'n':seg['U'],'H_over_I_upper':str(total/seg['U']),'decimal':dec(total/seg['U'])})
    # Decimal recomputation uses segment lengths instead of the per-d loop.
    with localcontext() as ctx:
        ctx.prec=85
        def D(q):return Decimal(q.numerator)/Decimal(q.denominator)
        for n,expected in [(best[1],best[0]),(cutoff,total/cutoff)]:
            ds=D(20*C[20])
            for seg in segments:
                count=max(0,min(n,seg['U'])-seg['L']+1)
                ds+=count*D(F(seg['upper']))
            assert abs(ds/n-D(expected))<Decimal('1e-75')
    tail_bound=F(tail['rational_upper']);constant=max(best[0],tail_bound,max(C.values()))
    safe=F(-(-(constant*10**6).numerator//(constant*10**6).denominator),10**6)
    assert constant<=safe<F('4.7475')
    result={'status':'certified consequence of source analytic lemmas and inherited small census',
        'all_dimension_constant_rational':str(constant),'all_dimension_constant_decimal':dec(constant),
        'safe_decimal_upper':dec(safe),'finite_chain_maximum_rational':str(best[0]),
        'finite_chain_maximum_decimal':dec(best[0]),'binding_dimension':best[1],
        'tail_upper_rational':str(tail_bound),'tail_upper_decimal':dec(tail_bound),'cutoff':cutoff,
        'base_dimension':20,'base_ratio':str(C[20]),'small_dimensions':{str(n):str(v) for n,v in C.items()},
        'small_block_choices':small_steps,'dimension_bounds_added':cutoff-20,'segments':segments,
        'checkpoints':checkpoints,'decimal_check_precision':85,
        'seconds':time.monotonic()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'nice':os.getpriority(os.PRIO_PROCESS,0),
        'certificate_sha256':{p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ('finite_certificate.json','blocks_certificate.json','tail_certificate.json','source_inputs/small_maxima.json')}}
    Path('ledger_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print('all dimension bound',dec(constant),'safe',dec(safe),'binding',best[1],flush=True)
    print('tail',dec(tail_bound),'cutoff',cutoff,'seconds',result['seconds'],flush=True)

if __name__=='__main__':main()
