import os,json,time,resource,importlib.util,sys
from pathlib import Path
root=Path.cwd();work=root/'computations/restriction_extension';work.mkdir(exist_ok=True)
sys.dont_write_bytecode=True
os.chdir(root/'computations/restriction_v4')
spec=importlib.util.spec_from_file_location('cert','certify.py');c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
row=json.loads(Path('proposals_improved_shift.json').read_text())[-1]
# Extend the finite range cheaply with the existing dual and one positive
# high-t multiplier. Its noise cost is small, while it controls x near zero.
row.update(L=8388609,U=33554432)
row['noise'].append(['8388608','4'])
row['alpha_padding']='1/1000'
os.chdir(work)
Path('proposal.json').write_text(json.dumps(row,indent=2)+'\n')
start=time.monotonic();cert=c.block(row)
Path('block_certificate.json').write_text(json.dumps(cert,indent=2)+'\n')
a=c.Q(1,45);n=row['U'];bound,terms=c.analytic(n,a)
tail={'n0':n,'a':str(a),'rational_upper':str(c.ceilq(bound)),'bound_interval':c.enc(bound),'terms_nats':[c.enc(z) for z in terms],'monotonicity':'Proposition 4.5, fixed a=1/45, n>=33554432'}
Path('tail_certificate.json').write_text(json.dumps(tail,indent=2)+'\n')
metrics={'seconds':time.monotonic()-start,'peak_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'nice':os.getpriority(os.PRIO_PROCESS,0)}
Path('metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
print('tail',float(c.ceilq(bound)),metrics,flush=True)
