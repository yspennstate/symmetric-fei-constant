"""Outward-interval checker. Floating discovery supplies no proof authority.

All input multipliers are rounded to fixed decimal rationals. Exact dimensions
check every spectral index. Dimension blocks check whole rational cells using
monotonicity, adaptively subdividing until every cell passes.
"""
import json, math, time, resource, os
from fractions import Fraction
from pathlib import Path
import mpmath as mp
iv=mp.iv
iv.dps=45
def Q(a,b=1):return iv.mpf(a)/b
def R(q):return Q(q.numerator,q.denominator)
def dec(x):return f'{x:.12f}'
def endpoint(v):return float(v.b)
def xlog(x):
    if x==0:return Q(0)
    assert x.a>0
    return x*x*iv.ln(x)
def hk(q):
    if q==1:return Q(0),Q(0)
    if q==Fraction(1,2):return 4*iv.ln(2),2*iv.ln(2)
    x=R(q);t=(1-x)/x;a=iv.sqrt(t)
    h=xlog(a+1)+xlog(abs(a-1))-t*iv.ln(t)
    ylog=lambda y: y*iv.ln(abs(y))
    k=a*(ylog(a+1)+ylog(a-1)-2*a*iv.ln(a))
    return h,k
def phi(v):
    term=Q(1);total=term;z=v*v/4
    j=0
    while True:
        j+=1;term*=z/(j*j);total+=term
        ratio=z/((j+2)**2)
        if ratio.b<Q(1).a:
            rem=term*z/((j+1)**2)/(1-ratio)
            if rem.b<(total*Q(1,10**55)).a:break
        assert j<10000
    lo=iv.exp(-v)*total;hi=iv.exp(-v)*(total+rem)
    return iv.mpf([lo.a,hi.b])
def params(row,lo):
    u=iv.mpf(dec(row['u']));B=iv.mpf(dec(row['mean']))
    noise=[(iv.mpf(dec(t)),iv.mpf(dec(w))) for t,w in row['noise'] if w>1e-8]
    assert B.a>=0 and all(w.a>=0 and t.a>0 for t,w in noise)
    C=Q(lo)-Q(1,6);z=C*u*u
    assert z.b<Q(1).a
    slope=1-2*iv.ln(u)+iv.ln(1-z)-2*z/(1-z)
    intercept=2*u/(1-z)
    cost=sum((w*phi(Q(lo-1)*(1-iv.exp(-2*t/lo))/4) for t,w in noise),Q(0))
    data=dict(u=dec(row['u']),B=dec(row['mean']),noise=[[dec(t),dec(w)] for t,w in row['noise'] if w>1e-8])
    return u,B,noise,slope,intercept,cost,data
def exact(row):
    n=row['n'];u,B,noise,slope,intercept,cost,data=params(row,n)
    values=[hk(Fraction(j,n)) for j in range(1,n+1)]
    requirements=[values[0][1]+slope-sum((w for t,w in noise),Q(0))]
    for s in range(1,n):
        x=Q(s,n)
        g=values[s-1][0]+values[s][1]-values[s-1][1]
        requirements.append(g+B*x-sum((w*iv.exp(-t*x) for t,w in noise),Q(0)))
    a=math.ceil(max(endpoint(v) for v in requirements)*10**9)+1
    alpha=Q(a,10**9)
    assert all(v.b<alpha.a for v in requirements)
    bound=(alpha+intercept-B*Q(n-1,2*n)+cost)/iv.ln(2)
    assert bound.b<Q(48677,10000).a,(n,str(bound))
    return dict(n=n,alpha=f'{a/10**9:.9f}',bound_bits=str(bound),bound_upper_float=endpoint(bound),coefficients_checked=len(requirements),**data)
def block(row):
    lo,hi=row['lo'],row['hi'];u,B,noise,slope,intercept,cost,data=params(row,lo)
    # Choose a rational constant with ample room below the target. Discovery
    # values only choose the candidate; every use is checked independently.
    pad=(4.8677-row['bound_bits'])*math.log(2)/3
    astr=f"{row['constant']+pad:.10f}";alpha=iv.mpf(astr)
    zero=hk(Fraction(1,lo))[1]+slope
    assert zero.b<(alpha+sum((w for t,w in noise),Q(0))).a
    M=(6+4*iv.ln(lo))/lo
    bound=(alpha+intercept-B*Q(lo-1,2*lo)+cost)/iv.ln(2)
    assert bound.b<Q(48677,10000).a
    cells=[];stack=[(Fraction(1,hi),Fraction(1,2)),(Fraction(1,2),Fraction(1))]
    count=0;minslack=1e9
    while stack:
        l,r=stack.pop();count+=1
        right=R(r)
        gap=alpha-B*right+sum((w*iv.exp(-t*right) for t,w in noise),Q(0))-hk(l)[0]-M
        if gap.a>0:
            cells.append([str(l),str(r)])
            minslack=min(minslack,float(gap.a))
        else:
            assert r-l>Fraction(1,10**12),(lo,hi,str(l),str(r),str(gap))
            mid=(l+r)/2;stack.extend([(l,mid),(mid,r)])
    Path(f'cells_{lo}_{hi}.json').write_text(json.dumps(cells)+'\n')
    return dict(lo=lo,hi=hi,alpha=astr,bound_bits=str(bound),bound_upper_float=endpoint(bound),cells=len(cells),evaluations=count,min_cell_slack_lower_float=minslack,**data)
def main():
    start=time.monotonic();out=dict(exact_dimensions=[],dimension_blocks=[])
    rows=json.loads(Path('duals.json').read_text())
    assert [r['n'] for r in rows]==list(range(16,256))
    for row in rows:
        result=exact(row);out['exact_dimensions'].append(result)
        if row['n']%20==0:print('certified n',row['n'],result['bound_upper_float'],flush=True)
        Path('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    blockranges=[(256,320),(320,500),(500,1000),(1000,2500),(2500,9999)]
    for lo,hi in blockranges:
        row=json.loads(Path(f'block_{lo}_{hi}.json').read_text())
        result=block(row);out['dimension_blocks'].append(result);print('certified block',result,flush=True)
        Path('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    out.update(status='certified',target_bits='4.8677',n0=16,precision_decimal_digits=iv.dps,seconds=time.monotonic()-start,peak_resident_bytes_macos=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,nice=os.getpriority(os.PRIO_PROCESS,0))
    Path('certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print('complete',out['seconds'],flush=True)
if __name__=='__main__':main()
