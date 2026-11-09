"""Second arithmetic implementation, with 90-digit outward intervals.

Replays every cell using squared-argument entropy formulas and an independent
primitive for the derivative envelope. Recomputes all finite coefficients.
"""
import gzip,json,os,resource,time
from fractions import Fraction as Q
from pathlib import Path
from mpmath import iv
iv.dps=90

def I(q):
    q=Q(q);return iv.mpf(q.numerator)/q.denominator
def ends(x):
    out=[]
    for sign,m,e,b in x._mpi_:
        z=Q(int(m))*Q(2)**int(e);out.append(-z if sign else z)
    return out
def lo(x):return ends(x)[0]
def hi(x):return ends(x)[1]
def mul_log(z):
    if lo(z)==hi(z)==0:return I(0)
    return z*iv.ln(z)
def pair(q):
    q=Q(q)
    if q==1:return I(0),I(0)
    if q==Q(1,2):return 4*iv.ln(2),2*iv.ln(2)
    x=I(q);t=(1-x)/x;p=iv.sqrt(t)
    a=p+1;b=p-1
    psi=(mul_log(a*a)+mul_log(b*b))/2-mul_log(t)
    k=p*(a*iv.ln(a*a)+b*iv.ln(b*b))/2-mul_log(t)
    return psi,k
def primitive(x):
    x=Q(x);v=1-2*x
    term=I(0) if v==0 else I(v)*iv.ln(abs(I(v)))
    return 2*I(x)+term-I(v)
def shift(l,r,L):
    h=Q(1,L);peak=Q(L-1,2*L)
    if peak<l:start=l
    elif peak>r:start=r
    else:start=peak
    return primitive(start+h)-primitive(start)
def phi(v):
    if hi(v)>100:return I(Q(24,23))/iv.sqrt(2*iv.pi*v)+iv.exp(-2*I(Q(23,48))**2*v)
    z=v*v/4;t=I(1);s=t
    for j in range(1,301):t=t*z/(j*j);s+=t
    remainder=t*z/(301*301)/(1-z/(302*302))
    assert hi(z/(302*302))<1
    return iv.exp(-v)*iv.mpf([s.a,(s+remainder).b])
def params(row,n):
    u=I(row['u']);beta=I(row['beta']);nw=[(I(t),I(w)) for t,w in row['noise']]
    assert lo(u)>0 and lo(beta)>=0 and all(lo(t)>0 and lo(w)>=0 for t,w in nw)
    v=1-(I(n)-I(Q(1,6)))*u*u;assert lo(v)>0
    slope=1-2*iv.ln(u)+iv.ln(v)-2*(1-v)/v;intercept=2*u/v
    cost=sum((w*phi(I(n-1)*(1-iv.exp(-2*t/n))/4) for t,w in nw),I(0))
    return beta,nw,slope,intercept,cost

def main():
    start=time.monotonic();finite=json.loads(Path('finite_certificate.json').read_text())['dimensions']
    source=json.loads(Path('source_inputs/finite_duals.json').read_text(),parse_float=str)
    checks=0
    for src,cert in zip(source,finite):
        n=src['n'];assert cert['n']==n
        row={'u':src['u'],'beta':src['mean'],'noise':src['noise']}
        beta,nw,slope,intercept,cost=params(row,n);alpha=I(cert['alpha'])
        assert hi(pair(Q(1,n))[1]+slope)<lo(alpha+sum((w for t,w in nw),I(0)))
        vals=[pair(Q(j,n)) for j in range(1,n+1)]
        for j in range(1,n):
            left=vals[j-1][0]+vals[j][1]-vals[j-1][1]
            right=alpha-beta*I(Q(j,n))+sum((w*iv.exp(-t*I(Q(j,n))) for t,w in nw),I(0))
            assert hi(left)<lo(right),(n,j)
        bound=(alpha+intercept-beta*I(Q(n-1,2*n))+cost)/iv.ln(2)
        assert hi(bound)<Q(cert['rational_upper']);checks+=n
    print('independent finite audit',checks,flush=True)
    rows=json.loads(Path('proposals_improved_shift.json').read_text())
    certs=json.loads(Path('blocks_certificate.json').read_text())['blocks'];cells=0;bounds=[]
    for row,cert in zip(rows,certs):
        n=row['L'];U=row['U'];assert (n,U)==(cert['L'],cert['U'])
        beta,nw,slope,intercept,cost=params(row,n);alpha=I(cert['alpha'])
        assert hi(pair(Q(1,n))[1]+slope)<lo(alpha+sum((w for t,w in nw),I(0)))
        with gzip.open(cert['covering'],'rt') as f:cover=json.load(f)
        cursor=Q(1,U)
        for a,b in cover:
            l=Q(a);r=Q(b);assert l==cursor and l<r<=1
            lhs=pair(l)[0]+shift(l,r,n)
            rhs=alpha-beta*I(r)+sum((w*iv.exp(-t*I(r)) for t,w in nw),I(0))
            assert hi(lhs)<lo(rhs),(n,U,l,r)
            cursor=r;cells+=1
        assert cursor==1 and len(cover)==cert['cells']
        bound=(alpha+intercept-beta*I(Q(n-1,2*n))+cost)/iv.ln(2)
        assert hi(bound)<Q(cert['rational_upper'])
        bounds.append({'L':n,'U':U,'bound_upper_rational':str(hi(bound))})
        print('independent block audit',n,U,len(cover),flush=True)
    # Point arithmetic with integral quadrature cross-checks the modified tail.
    # The proof enclosure itself remains the outward interval in tail_certificate.
    import mpmath as mp
    mp.mp.dps=85
    tail=json.loads(Path('tail_certificate.json').read_text());n=tail['n0'];N=mp.mpf(n-1);a=mp.mpf(1)/30;r=mp.mpf(23)/48
    c=N/n*(-mp.expm1(-2*a))/(2*a);vT=N*(-mp.expm1(-2*a))/4;v20=N*(-mp.expm1(-mp.mpf(40)/n))/4
    def integral_phi(v):return 2/mp.pi*mp.quad(lambda t:mp.exp(-2*v*mp.sin(t)**2),[0,mp.pi/2])
    b=mp.sqrt(2/(mp.pi*N));z=40*r*r;A=-mp.mpf(1)/2-mp.log(a)+a
    result=(b*(3+mp.log1p(1/b**2)+2*mp.log(2))+(6+4*mp.log(n))/n+mp.mpf(14)/5+mp.mpf(273)/(100*n)
        +48/(23*mp.sqrt(40*mp.pi))+mp.exp(-z)/z+integral_phi(20*c)*mp.log(1/c)
        +mp.exp(a)*A*(24/(23*mp.sqrt(2*mp.pi*vT))+mp.exp(-2*r*r*vT))
        +mp.mpf(89)/50*integral_phi(v20))/mp.log(2)
    interval=[mp.mpf(Q(x).numerator)/Q(x).denominator for x in tail['bound_interval']]
    assert interval[0]<=result<=interval[1]
    ans={'status':'passed','finite_coefficient_inequalities':checks,'block_cells':cells,'blocks':len(rows),
        'interval_digits':iv.dps,'tail_quadrature_decimal':str(result),'tail_quadrature_digits':85,
        'seconds':time.monotonic()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        'nice':os.getpriority(os.PRIO_PROCESS,0),'block_bounds':bounds}
    Path('audit.json').write_text(json.dumps(ans,indent=2)+'\n');print('independent audit passed',ans['seconds'],flush=True)

if __name__=='__main__':main()
