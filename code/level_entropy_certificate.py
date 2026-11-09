"""Outward-rounded interval certificate for Lemma 3.2: L <= (673/250) I for every Boolean f.

Run: python3 level_entropy_certificate.py   (writes level_entropy_certificate.json, about 3 s)
"""
import json,time
from pathlib import Path
import mpmath as mp

iv=mp.iv
iv.dps=45
one=iv.mpf(1);two=iv.mpf(2);ln2=iv.log(two)
r=iv.mpf(147)/100
C=iv.mpf(673)/250 # 2.692 exactly
A=iv.log(r)/ln2

def upper_on_interval(lo,hi):
    q=iv.mpf([lo,hi]);p=iv.exp(-q*ln2);e=4*p*(1-p)
    h=-e*iv.log(e)-(1-e)*iv.log(1-e)
    return A+(h-e*iv.log(r-1))/(2*p*q*ln2)

def prove():
    start=time.monotonic()
    large_p=A+(iv.log(two)-iv.log(r-1))/ln2
    assert large_p.b<C.a
    cutoff=32
    tail=A+2+2*(1-iv.log(iv.mpf(4))-iv.log(r-1)-iv.log(1-iv.mpf(2)**(-cutoff)))/(cutoff*ln2)
    assert tail.b<C.a
    stack=[(2.,32.)];leaves=[];max_upper=iv.mpf(0)
    while stack:
        lo,hi=stack.pop();b=upper_on_interval(lo,hi)
        if b.b<C.a:
            leaves.append([lo,hi,str(b.b)])
            if b.b>max_upper:max_upper=b.b
        else:
            assert hi-lo>2**(-30),'certificate failed to converge'
            mid=(lo+hi)/2
            stack.extend([(mid,hi),(lo,mid)])
    leaves.sort()
    assert leaves[0][0]==2 and leaves[-1][1]==32
    assert all(leaves[j][1]==leaves[j+1][0] for j in range(len(leaves)-1))
    result={'status':'certified','C_exact':'673/250','C_decimal':'2.692',
            'geometric_base_exact':'147/100','geometric_a_interval':str(A),
            'large_p_region':'1/4 <= p <= 1/2','large_p_upper_interval':str(large_p),
            'small_p_region':'0 < p <= 2^-32','small_p_upper_interval':str(tail),
            'compact_region_q':'2 <= -log2(p) <= 32',
            'maximum_leaf_upper_interval':str(max_upper),
            'interval_precision_decimal_digits':iv.dps,'leaf_count':len(leaves),
            'elapsed_seconds':time.monotonic()-start,'intervals':leaves}
    Path('level_entropy_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='intervals'},indent=2))

if __name__=='__main__':prove()
