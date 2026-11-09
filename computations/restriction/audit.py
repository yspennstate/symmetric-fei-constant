"""Independent Decimal and exact-integer audits; not substituted for proof checks."""
from decimal import Decimal as D, localcontext
from pathlib import Path
from fractions import Fraction
import math,json,time,resource,hashlib

def hk_dec(x):
    if x==1:return D(0),D(0)
    if x==D('.5'):return 4*D(2).ln(),2*D(2).ln()
    t=(1-x)/x;a=t.sqrt()
    h=(a+1)**2*(a+1).ln()+(a-1)**2*abs(a-1).ln()-t*t.ln()
    k=a*((a+1)*(a+1).ln()+(a-1)*abs(a-1).ln()-2*a*a.ln())
    return h,k
def phi_dec(v):
    term=total=D(1);j=0
    while True:
        j+=1;term*=v*v/(4*j*j);total+=term
        if j>v and term<total*D('1e-80'):break
    return (-v).exp()*total
def decimal_audit(cert):
    points=0;worst_gap=D('Infinity');max_error=D(0);bounds=[]
    for row in cert['exact_dimensions']:
        n=row['n'];u=D(row['u']);B=D(row['B']);alpha=D(row['alpha']);C=D(n)-D(1)/6
        pairs=[hk_dec(D(s)/n) for s in range(1,n+1)]
        y=C*u*u;slope=1-2*u.ln()+(1-y).ln()-2*y/(1-y);intercept=2*u/(1-y)
        noise=[(D(t),D(w)) for t,w in row['noise']]
        gap=alpha+sum(w for t,w in noise)-pairs[0][1]-slope
        assert gap>0
        worst_gap=min(worst_gap,gap);points+=1
        for s in range(1,n):
            x=D(s)/n;g=pairs[s-1][0]+pairs[s][1]-pairs[s-1][1]
            gap=alpha-B*x+sum(w*(-t*x).exp() for t,w in noise)-g
            assert gap>0
            worst_gap=min(worst_gap,gap);points+=1
        bound=(alpha+intercept-B*D(n-1)/(2*n)+sum(w*phi_dec(D(n-1)*(1-(-2*t/n).exp())/4) for t,w in noise))/D(2).ln()
        lo,hi=[D(s.strip()) for s in row['bound_bits'][1:-1].split(',')]
        assert lo<bound<hi
        bounds.append((bound,n));max_error=max(max_error,abs(bound-(lo+hi)/2))
    block_cells=0
    for row in cert['dimension_blocks']:
        L,U=row['lo'],row['hi'];alpha=D(row['alpha']);B=D(row['B']);M=(6+4*D(L).ln())/L
        noise=[(D(t),D(w)) for t,w in row['noise']]
        cells=[(Fraction(l),Fraction(r)) for l,r in json.loads(Path(f'cells_{L}_{U}.json').read_text())]
        ordered=sorted(cells)
        assert ordered[0][0]==Fraction(1,U) and ordered[-1][1]==1
        assert all(a[1]==b[0] for a,b in zip(ordered,ordered[1:]))
        for l,r in cells:
            dl=D(l.numerator)/l.denominator;dr=D(r.numerator)/r.denominator
            gap=alpha-B*dr+sum(w*(-t*dr).exp() for t,w in noise)-hk_dec(dl)[0]-M
            assert gap>0;block_cells+=1
    return dict(decimal_precision_digits=75,coefficients_checked=points,block_cells_checked=block_cells,
        min_coefficient_slack_nats=str(worst_gap),max_midpoint_difference_bits=str(max_error),
        largest_exact_dimension_bound_bits=str(max(bounds)[0]),largest_exact_dimension_bound_dimension=max(bounds)[1])

def integer_walsh_audit():
    import small_dimensions as sm
    counts=full_walsh=0;max_tangent_error=0.;max_zero_to_tangent=-1e100
    def xl(x):return x*x*math.log(abs(x)) if x else 0.
    def pair(a,b):return xl(a+b)+xl(a-b)-2*xl(a)-2*xl(b)
    for n in range(2,11):
        N=n-1;mat=sm.kraw(n);restmat=sm.kraw(N);cm=[math.comb(N,s) for s in range(n)]
        for mask in range(1,1<<n):
            f=[1]+[-1 if mask>>(w-1)&1 else 1 for w in range(1,n+1)]
            nums=sm.transform(f,mat)
            assert sum(math.comb(n,s)*v*v for s,v in enumerate(nums))==4**n
            weights=[cm[s]*nums[s+1]**2 for s in range(n)];T=sum(weights)
            assert 2*sum(s*w for s,w in enumerate(weights))>=N*T
            edge=sum(cm[s] for s in range(n) if f[s]!=f[s+1])
            assert T==edge*(1<<(n+1))
            assert nums[0]**2+n*nums[1]**2<=4**n
            assert abs(nums[1])*(1<<n)<=T
            mu=[v/T for v in weights];a=[v/(1<<n) for v in nums];I=T/4**n
            L0=pair(a[0],a[1])/I
            u=.4/math.sqrt(n);C=n-1/6;y=C*u*u
            slope=1-2*math.log(u)+math.log1p(-y)-2*y/(1-y);intercept=2*u/(1-y)
            diff=L0-(intercept+slope*mu[0]);max_zero_to_tangent=max(max_zero_to_tangent,diff)
            assert diff<1e-11
            def hkf(x):
                if x==1:return 0.,0.
                if x==.5:return 4*math.log(2),2*math.log(2)
                t=(1-x)/x;r=math.sqrt(t)
                yl=lambda z:z*math.log(abs(z)) if z else 0.
                return pair(r,1),r*(yl(r+1)+yl(r-1)-2*yl(r))
            hk=[hkf(s/n) for s in range(1,n+1)]
            delta=sum(cm[s]*pair(a[s],a[s+1]) for s in range(n))/I
            tangent=L0+mu[0]*hk[0][1]+sum(mu[s]*(hk[s-1][0]+hk[s][1]-hk[s-1][1]) for s in range(1,n))
            assert delta<=tangent+1e-10
            # Independent direct restriction-entropy calculation.
            pp=sm.transform(f[:-1],restmat);mm=sm.transform(f[1:],restmat)
            ent=lambda nums,d:-2*sum(math.comb(d,s)*xl(v/(1<<d)) for s,v in enumerate(nums))
            direct=(ent(nums,n)-(ent(pp,N)+ent(mm,N))/2)/I
            max_tangent_error=max(max_tangent_error,abs(delta-direct))
            assert abs(delta-direct)<1e-10
            if n<=7:
                literal=[sum(f[x.bit_count()] *(-1)**((x&((1<<s)-1)).bit_count()) for x in range(1<<n)) for s in range(n+1)]
                assert literal==nums;full_walsh+=1
            counts+=1
    return dict(nonconstant_tables_modulo_output_negation=counts,full_cube_walsh_checks=full_walsh,
        maximum_pair_vs_direct_entropy_error_nats=max_tangent_error,maximum_zero_level_minus_tangent_nats=max_zero_to_tangent)

def assembly(cert):
    from mpmath import iv
    iv.dps=70
    K=12-iv.mpf(9)*iv.ln(3)/(2*iv.ln(2))
    def parse(text):
        a,b=text[1:-1].split(',')
        return iv.mpf([iv.mpf(a.strip()).a,iv.mpf(b.strip()).b])
    assert all(parse(r['bound_bits']).b<K.a for r in cert['exact_dimensions']+cert['dimension_blocks'])
    large=json.loads(Path('cutoff_certificate.json').read_text())
    assert parse(large['bound_bits_per_unit_coordinate_influence']).b<K.a
    small=json.loads(Path('small_certificate.json').read_text())
    assert [r['n'] for r in small['rows']]==list(range(1,16))
    for row in small['rows']:
        if row['n']==3:
            signature={int(p):Fraction(v) for p,v in row['exact_ratio_prime_log_coefficients'].items()}
            assert signature=={2:Fraction(12),3:Fraction(-9,2)},signature
        else:assert iv.mpf(row['ratio_decimal_enclosure'][1]).b<K.a
    # Integer table counts and complete interval coverage are checked separately.
    assert [r['n'] for r in cert['exact_dimensions']]==list(range(16,256))
    assert [(r['lo'],r['hi']) for r in cert['dimension_blocks']]==[(256,320),(320,500),(500,1000),(1000,2500),(2500,9999)]
    return dict(status='all dimensions covered; sharp restriction constant certified subject to analytic proof review',
        analytic_cutoff=16,combined_cutoff=1,sharp_restriction_constant_bits=str(K),
        small_sign_tables=sum(r['tables'] for r in small['rows']),large_dimension_bound_bits=large['bound_bits_per_unit_coordinate_influence'])

def main():
    start=time.monotonic();cert=json.loads(Path('certificate.json').read_text());assert cert['status']=='certified'
    with localcontext() as ctx:
        ctx.prec=75
        out=dict(decimal_audit=decimal_audit(cert))
    out['integer_walsh_audit']=integer_walsh_audit();out['assembly']=assembly(cert)
    out['status']='passed';out['seconds']=time.monotonic()-start
    out['peak_resident_bytes_macos']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    Path('audit.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
