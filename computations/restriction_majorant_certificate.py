"""Outward interval certificate for a scalar majorant; no scipy needed."""
import mpmath as mp
import math,json
from pathlib import Path
mp.iv.dps=45
iv=mp.iv

def pt(x):return iv.mpf(str(x))

def F(x):
    if x==.5:return 4*iv.ln(2)
    xx=pt(x);t=(1-xx)/xx;a=iv.sqrt(t)
    return (a+1)**2*iv.ln(a+1)+(a-1)**2*iv.ln(abs(a-1))-t*iv.ln(t)

def e1_lower(z):
    if float(z.b)<=1:
        # E1's alternating series, stopped after its eighth (negative) term.
        sm=-pt('.578')-iv.ln(z)
        for k in range(1,9):sm+=(-1)**(k+1)*z**k/(k*math.factorial(k))
        return sm
    # Three-dimensional polynomial variational lower bound for the Stieltjes transform.
    return iv.exp(-z)*(z*z+8*z+11)/(z**3+9*z*z+18*z+6)

def U_lower(x):
    xx=pt(x)
    return e1_lower(40*xx)+pt('5.53')-pt('5.46')*xx+pt('1.78')*iv.exp(-20*xx)

def main():
    # Gamma upper bound from the Euler--Maclaurin remainder, using m=10.
    gamma_upper=sum((iv.mpf(1)/k for k in range(1,11)),iv.mpf(0))-iv.ln(10)-iv.mpf(1)/20+iv.mpf(1)/1200
    assert float(gamma_upper.b)<.578
    # For x<=1e-4: F<=-ln x+3, E1(40x)>=-gamma-ln(40x), e^-20x>=1-20x.
    small_slack=pt('5.53')+pt('1.78')-pt('.578')-iv.ln(40)-3-(pt('5.46')+20*pt('1.78'))*pt('.0001')
    assert float(small_slack.a)>0
    nodes=[.0001]+[math.exp(math.log(.0001)+j*math.log(100)/2000) for j in range(1,2000)]+[.01]
    nodes += [.01+k/10000 for k in range(1,9901)]
    nodes[-1]=1.
    worst=1e9;where=None
    for idx,(l,r) in enumerate(zip(nodes,nodes[1:])):
        # True U is decreasing, and F is decreasing. U_lower(r)>=F(l) proves the whole cell.
        slack=U_lower(r)-F(l)
        low=float(slack.a)
        assert low>0,(l,r,str(slack))
        if low<worst:worst=low;where=[l,r]
        if idx%3000==0:print('cells_checked',idx,'running_minimum_slack',worst,flush=True)
    # E E1(40X) <= integral_40^infinity exp(-t/2)I0(t/2)dt/t.
    # Split the Bessel angular integral at theta=1/2. On [0,1/2],
    # sin(theta)>=2 sin(1/2)*theta >= (23/24)*theta.
    # On [1/2,pi/2], sin(theta)>=sin(1/2)>=23/48.
    # Integrate the resulting Gaussian term and E1 upper bound.
    sin_lower=iv.mpf(23)/48
    zz=40*sin_lower**2
    tail=1/(iv.sqrt(40*iv.pi)*sin_lower)+iv.exp(-zz)/zz
    # At t=20, exact positive series I0(10), e^-10; bound with a geometric tail.
    term=iv.mpf(1);sm=term
    for k in range(1,61):
        term*=iv.mpf(25)/(k*k);sm+=term
    next_term=term*iv.mpf(25)/(61*61)
    ratio=iv.mpf(25)/(62*62)
    upperI=sm+next_term/(1-ratio)
    B=iv.exp(-10)*upperI
    bound=(tail+pt('5.53')-pt('5.46')/2+pt('1.78')*B)/iv.ln(2)
    out=dict(cells=len(nodes)-1,min_cell_slack=math.nextafter(worst,-math.inf),worst_cell=where,
             small_endpoint_slack_lower=math.nextafter(float(small_slack.a),-math.inf),gamma_upper=math.nextafter(float(gamma_upper.b),math.inf),
             bessel10_upper=math.nextafter(float(B.b),math.inf),tail_upper=math.nextafter(float(tail.b),math.inf),bound_bits_upper=math.nextafter(float(bound.b),math.inf),
             pointwise_majorant='E1(40x)+5.53-5.46x+1.78 exp(-20x)',
             status='outward-interval certificate; continuum moment constraints required')
    assert float(bound.b)<4.637
    Path('majorant_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)

if __name__=='__main__':main()
