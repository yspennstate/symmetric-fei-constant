"""Fresh scalar certificate: exact rational mesh, Taylor E1 and Jensen tail.

Does not import or run any w18 code. mpmath intervals round outwards.
For z<=4 the integral Taylor remainder gives an E1 lower bound after
40 terms. For z>4 Jensen gives E1(z)>=exp(-z)/(z+1).
"""
import json, math, time
from pathlib import Path
import mpmath as mp

iv=mp.iv
iv.dps=55
def rat(p,q=1): return iv.mpf(p)/q
def lower(x): return float(x.a)
def upper(x): return float(x.b)
def xlog(x): return x*x*iv.ln(abs(x))
def potential(x):
    # Only rational point intervals, no cell straddles a singularity.
    if x == rat(1,2): return 4*iv.ln(2)
    if x == rat(1): return rat(0)
    t=(1-x)/x; a=iv.sqrt(t)
    return xlog(a+1)+xlog(a-1)-t*iv.ln(t)

GAMMA=rat(57722,100000)
def e1_lower(z):
    if upper(z)<=4:
        ans=-GAMMA-iv.ln(z)
        for j in range(1,41):
            ans+=(-1)**(j+1)*z**j/(j*math.factorial(j))
        return ans
    return iv.exp(-z)/(z+1)
def majorant_lower(x):
    return e1_lower(40*x)+rat(553,100)-rat(273,50)*x+rat(89,50)*iv.exp(-20*x)

def main():
    start=time.time()
    gamma_upper=sum((rat(1,k) for k in range(1,21)),rat(0))-iv.ln(20)-rat(1,40)+rat(1,4800)
    assert upper(gamma_upper)<lower(GAMMA)
    # Entire (0,1e-4] covered analytically.
    endpoint=rat(553,100)+rat(89,50)-GAMMA-iv.ln(40)-3-(rat(273,50)+rat(178,5))*rat(1,10000)
    assert lower(endpoint)>0
    # Exact decimal-rational mesh; small region uses 1e-6, rest 1e-4.
    nodes=[rat(j,1000000) for j in range(100,10001)]
    nodes += [rat(j,10000) for j in range(101,10001)]
    worst=float('inf'); worst_cell=None
    for i,(left,right) in enumerate(zip(nodes,nodes[1:])):
        gap=majorant_lower(right)-potential(left)
        v=lower(gap)
        assert v>0,(i,str(left),str(right),str(gap))
        if v<worst:worst=v;worst_cell=[str(left),str(right)]
        if i%4000==0:print('certified_cells',i,'min_slack',worst,flush=True)
    # Independent final arithmetic: deliberately overbound B(10) by a rational.
    term=rat(1); sm=term
    for j in range(1,81):term*=rat(25,j*j);sm+=term
    bupper=iv.exp(-10)*(sm+term*rat(25,81**2)/(1-rat(25,82**2)))
    assert upper(bupper)<0.127833337164
    z=40*rat(23,48)**2
    J=rat(48,23)/iv.sqrt(40*iv.pi)+iv.exp(-z)/z
    final=(J+rat(14,5)+rat(89,50)*bupper)/iv.ln(2)
    assert upper(final)<4.636423666
    out=dict(status='proved by outward intervals and analytic endpoint bounds',
        cells=len(nodes)-1,min_slack_lower=worst,worst_cell=worst_cell,
        small_endpoint_slack_lower=lower(endpoint),gamma_upper=upper(gamma_upper),
        bessel10_upper=upper(bupper),J_upper=upper(J),constant_bits_upper=upper(final),
        seconds=time.time()-start)
    Path('independent_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':main()
