import json
from pathlib import Path
rows=json.loads(Path('certified_maxima.json').read_text())
lines=['// Generated from certified_maxima.json; verifier recomputes all values.','static const uint32_t REFERENCE[25][6] = {']
for n in range(25):
 vals=[0]*6
 for r in rows:
  if r['n']==n:
   t=r['maximizers_mod_output_complement'][0];vals[r['k']]=sum(int(v)<<w for w,v in enumerate(t))
 lines.append('{'+','.join(map(str,vals))+'},')
lines.append('};')
Path('references.h').write_text('\n'.join(lines)+'\n')
