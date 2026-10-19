#!/usr/bin/env python3
"""Execute the released fixed analyses into a new directory, without old outputs."""
import argparse, datetime, json, os, pathlib, subprocess, sys

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--out',type=pathlib.Path,required=True)
    p.add_argument('--scope',choices=['deterministic','all'],default='all')
    p.add_argument('--r-command',help='Pinned R executable or launcher (not Rscript)')
    a=p.parse_args()
    if a.scope=='all' and not a.r_command:p.error('--scope all requires --r-command; see sensitivity/README.md')
    root=pathlib.Path(__file__).resolve().parent
    out=a.out.resolve()
    if out.exists():p.error('Output directory must not already exist')
    out.mkdir(parents=True)
    env=os.environ.copy()
    env.update(PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    jobs=[('descriptive','descriptive/run_public.py','--out'),('cr2','cr2/run_public.py','--out'),('exploration','exploration/run_all.py','--out'),('fixed_model_reconstruction','fixed_model_reconstruction/run_reconstruction.py','--output-dir'),('support_selection','support_selection/run_selection.py','--out')]
    if a.scope=='all':jobs.append(('sensitivity','sensitivity/run_public.py','--out'))
    done=[]
    for name,script,flag in jobs:
        cmd=[sys.executable,'-B',str(root/script),flag,str(out/name)]
        if name=='sensitivity':cmd+=['--r-command',a.r_command]
        print('Running '+name,flush=True)
        with (out/(name+'.log')).open('w',encoding='utf-8') as log:
            subprocess.run(cmd,cwd=root,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        done.append(name)
    receipt={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':a.scope,'fresh_execution_modules':done,'tost_executed':a.scope=='all','manuscript_comparison':False,'original_source_collection_reproduced':False}
    (out/'EXECUTION_RECEIPT.json').write_text(json.dumps(receipt,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(receipt,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
