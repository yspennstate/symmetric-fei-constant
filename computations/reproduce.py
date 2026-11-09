"""Rebuild and audit the dimension-24 seed, continuous extension and ledger."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[1]
os.chdir(root)
for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
if os.getpriority(os.PRIO_PROCESS, 0) < 10:
    os.setpriority(os.PRIO_PROCESS, 0, 10)

def run(args, cwd=root):
    print('Running', ' '.join(map(str, args)), flush=True)
    subprocess.run(args, cwd=cwd, check=True)

census = root / 'computations/block_census_24'
with tempfile.TemporaryDirectory(dir=root, prefix='.certificate-') as scratch:
    binary = Path(scratch) / 'enumerate'
    run(['c++', '-O3', '-std=c++17', census / 'enumerate.cpp', '-o', binary])
    run([binary, '24'], census)
    run([sys.executable, 'audit.py'], census)
    run([sys.executable, 'prepare_check.py'], census)
    binary = Path(scratch) / 'check_census'
    run(['c++', '-O3', '-std=c++17', census / 'check_census.cpp', '-o', binary])
    run([binary], census)
for script in ('certify_extension.py', 'check_extension.py', 'rebuild_ledger.py', 'audit_ledger.py'):
    run([sys.executable, root / 'computations' / script])
print('All certificates and audits passed.', flush=True)
