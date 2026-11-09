"""Outward-interval certification of finite duals, continuous blocks and tail.

The LP is not imported. Decimal proposals are treated as exact rationals.
Only Proposition 4.4/4.5 and the derivative lemma in paper/main.tex are used.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS'):
    os.environ[name]='1'
import argparse, gzip, hashlib, importlib.util, json, resource, sys, time
from fractions import Fraction as Q
from pathlib import Path
sys.dont_write_bytecode=True
spec=importlib.util.spec_from_file_location('reference',Path('source_inputs/reference_certificate.py'))
ref=importlib.util.module_from_spec(spec);spec.loader.exec_module(ref)
iv=ref.iv;iv.dps=65
I=ref.Q;R=ref.R

def endpoint(raw):
    sign,man,exp,_=raw;q=Q(int(man))*Q(2)**int(exp)
    return -q if sign else q
def low(x):return endpoint(x._mpi_[0])
def high(x):return endpoint(x._mpi_[1])
def ceilq(x,places=15):
    a=high(x)*10**places
    return Q(-(-a.numerator//a.denominator),10**places)
def enc(x):return [str(low(x)),str(high(x))]
def phi(v):
    if high(v)<=100:return ref.phi(v)
    return I(24,23)/iv.sqrt(2*iv.pi*v)+iv.exp(-2*I(23,48)**2*v)

def setup(row,L):
    u=R(Q(row['u']));beta=R(Q(row['beta']))
    nw=[(R(Q(t)),R(Q(w))) for t,w in row['noise']]
    assert low(u)>0 and low(beta)>=0
    assert all(low(t)>0 and low(w)>=0 for t,w in nw)
    C=I(L)-I(1,6);z=C*u*u;assert high(z)<1
    j=1-2*iv.ln(u)+iv.ln(1-z)-2*z/(1-z)
    d=2*u/(1-z)
    cost=sum((w*phi(I(L-1)*(1-iv.exp(-2*t/L))/4) for t,w in nw),I(0))
    return beta,nw,j,d,cost

def finite():
    start=time.monotonic();out=[]
    for old in json.loads(Path('source_inputs/finite_duals.json').read_text(),parse_float=str):
        n=old['n'];row={'u':old['u'],'beta':old['mean'],'noise':old['noise']}
        beta,nw,j,d,cost=setup(row,n)
        values=[ref.hk(Q(s,n)) for s in range(1,n+1)]
        req=[values[0][1]+j-sum((w for t,w in nw),I(0))]
        for s in range(1,n):
            x=I(s,n)
            req.append(values[s-1][0]+values[s][1]-values[s-1][1]+beta*x-sum((w*iv.exp(-t*x) for t,w in nw),I(0)))
        maxreq=max(high(z) for z in req)
        alpha=Q(-(-(maxreq*10**35).numerator//(maxreq*10**35).denominator),10**35)+Q(1,10**33)
        assert all(high(v)<low(R(alpha)) for v in req)
        bound=(R(alpha)+d-beta*I(n-1,2*n)+cost)/iv.ln(2)
        out.append({'n':n,'alpha':str(alpha),'rational_upper':str(ceilq(bound)),
            'bound_interval':enc(bound),'coefficients_checked':n})
    assert [r['n'] for r in out]==list(range(16,256))
    result={'dimensions':out,'coefficients_checked':sum(r['coefficients_checked'] for r in out),'seconds':time.monotonic()-start,'digits':iv.dps}
    Path('finite_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print('finite certified',result['coefficients_checked'],'inequalities',result['seconds'],'seconds',flush=True)

def primitive(z):
    if z==0:return I(0)
    v=R(z);return 4*v-2*v*iv.ln(2*abs(v))

def shift_envelope(l,r,L):
    # B(y)=2-2 log|1-2y| is nonnegative on [0,1+1/L], L>=256.
    # The integral over a fixed length is maximal when centred at 1/2.
    # Clamp that maximizer to [l,r] to cover every possible cell starting point.
    h=Q(1,L);x=min(r,max(l,Q(1,2)-h/2))
    return primitive(x+h-Q(1,2))-primitive(x-Q(1,2))

def block(row,replay=False):
    start=time.monotonic();L=row['L'];U=row['U'];assert 256<=L<=U
    beta,nw,j,d,cost=setup(row,L)
    aq=Q(row['alpha'])+Q(row['alpha_padding']);alpha=R(aq)
    zero=alpha+sum((w for t,w in nw),I(0))-ref.hk(Q(1,L))[1]-j
    assert low(zero)>0,(L,U,'zero atom')
    local=row.get('shift_bound')=='improved_integrated_derivative_envelope'
    def gap(l,r):
        M=shift_envelope(l,r,L) if local else (6+4*iv.ln(L))/L
        return alpha-beta*R(r)+sum((w*iv.exp(-t*R(r)) for t,w in nw),I(0))-ref.hk(l)[0]-M
    name=f'cells_{L}_{U}.json.gz';cells=[];minslack=None
    if replay:
        with gzip.open(name,'rt') as f:data=json.load(f)
        cursor=Q(1,U)
        for a,b in data:
            l=Q(a);r=Q(b);assert l==cursor and l<r<=1
            z=gap(l,r);assert low(z)>0,(L,U,l,r,'replay')
            minslack=min(minslack,low(z)) if minslack is not None else low(z)
            cells.append([a,b]);cursor=r
        assert cursor==1
    else:
        stack=[(Q(1,U),Q(1,2),0),(Q(1,2),Q(1),0)]
        while stack:
            l,r,depth=stack.pop();z=gap(l,r)
            if low(z)>0:
                cells.append([str(l),str(r)])
                minslack=min(minslack,low(z)) if minslack is not None else low(z)
            else:
                assert depth<55,(L,U,l,r,'unresolved cell',str(z))
                m=(l+r)/2;stack.extend([(l,m,depth+1),(m,r,depth+1)])
        cells.sort(key=lambda a:Q(a[0]));cursor=Q(1,U)
        for a,b in cells:assert Q(a)==cursor;cursor=Q(b)
        assert cursor==1
        with gzip.open(name,'wt') as f:json.dump(cells,f,separators=(',',':'))
    bound=(alpha+d-beta*I(L-1,2*L)+cost)/iv.ln(2)
    result={'L':L,'U':U,'alpha':str(aq),'rational_upper':str(ceilq(bound)),
        'bound_interval':enc(bound),'cells':len(cells),'min_slack_nats':str(minslack),
        'zero_slack_nats':str(low(zero)),'covering':name,'shift_bound':row.get('shift_bound','global'),
        'seconds':time.monotonic()-start}
    print('block certified',L,U,float(ceilq(bound)),len(cells),'cells',result['seconds'],'seconds',flush=True)
    return result

def analytic(n,a=Q(1,30)):
    # Proposition 4.5 with a rational 0<a<=1/5; requires an>=40.
    assert 0<a<=Q(1,5) and a*n>=40
    N=I(n-1);a=R(a);r=I(23,48)
    c=N/n*(1-iv.exp(-2*a))/(2*a)
    vT=N*(1-iv.exp(-2*a))/4;v20=N*(1-iv.exp(-I(40,n)))/4
    b=iv.sqrt(2/(iv.pi*N));z=40*r*r
    A=-I(1,2)-iv.ln(a)+a
    terms=[b*(3+iv.ln(1+1/(b*b))+2*iv.ln(2)),(6+4*iv.ln(n))/n,
        I(14,5),I(273,100*n),I(48,23)/iv.sqrt(40*iv.pi)+iv.exp(-z)/z,
        ref.phi(20*c)*iv.ln(1/c),iv.exp(a)*A*(I(24,23)/iv.sqrt(2*iv.pi*vT)+iv.exp(-2*r*r*vT)),
        I(89,50)*ref.phi(v20)]
    return sum(terms,I(0))/iv.ln(2),terms

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--replay',action='store_true');parser.add_argument('--finite-only',action='store_true');args=parser.parse_args()
    t=time.monotonic()
    finite()
    if args.finite_only:return
    rows=json.loads(Path('proposals_improved_shift.json').read_text());out=[]
    for row in rows:
        out.append(block(row,args.replay))
        Path('blocks_replay.json' if args.replay else 'blocks_certificate.json').write_text(json.dumps({'blocks':out,'complete':len(out)==len(rows),'digits':iv.dps},indent=2)+'\n')
    cutoff=rows[-1]['U'];bound,terms=analytic(cutoff)
    result={'n0':cutoff,'a':'1/30','rational_upper':str(ceilq(bound)),'bound_interval':enc(bound),'terms_nats':[enc(z) for z in terms],
        'monotonicity':'same termwise argument as Proposition 4.5, fixed a=1/30; n>=1200'}
    Path('tail_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    Path('certification_metrics.json').write_text(json.dumps({'seconds':time.monotonic()-t,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'nice':os.getpriority(os.PRIO_PROCESS,0),'replay':args.replay},indent=2)+'\n')
    print('tail certified',float(ceilq(bound)),'from',cutoff,flush=True)

if __name__=='__main__':main()
