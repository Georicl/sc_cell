"""顺序重跑既有计算核验，保存本轮复核日志。"""
from pathlib import Path
import subprocess
import sys
import json

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'results/Wu2021_multislice_CAF_T_v1/previous_result_review'
OUT.mkdir(parents=True,exist_ok=True)
scripts=['audit_pipeline_to_qc.py','validate_known_composition.py','validate_reference_comparison.py','validate_independent_validation.py','validate_tcell_validation.py','validate_tcell_label_review.py','validate_coarse_tcells.py']
records=[]
for name in scripts:
    log=OUT/(Path(name).stem+'.log')
    with log.open('w') as f:
        result=subprocess.run([sys.executable,str(ROOT/'src/sc_sncell_github'/name)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    records.append(dict(script=name,exit_code=result.returncode,log=str(log.relative_to(ROOT))))
    print(name,'PASS' if result.returncode==0 else 'FAIL',flush=True)
(OUT/'review_status.json').write_text(json.dumps(records,indent=2)+'\n')
assert all(x['exit_code']==0 for x in records)
