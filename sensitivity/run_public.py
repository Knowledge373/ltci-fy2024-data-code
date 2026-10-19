#!/usr/bin/env python3
"""Fresh fixed-support computation: six GLMs and all 1,600 TOST boundaries by default."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import argparse,csv,datetime,hashlib,importlib.util,io,json,math,shutil,subprocess,sys,time
from pathlib import Path
import numpy as np
sys.dont_write_bytecode=True
ROOT=Path(__file__).resolve().parent
KEYS=['P1_k1','P1_k0p5','P2_k1','P2_k0p5','P3_k1','P3_k0p5']
SOURCE_SHA='ea25bf07ed1ff42626ca73719ca252af7f563ca35cd19029139863cdbca468f2'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def bh(b):return hashlib.sha256(b).hexdigest()
def dump(p,z):Path(p).write_text(json.dumps(z,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def rows(p):
 with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def load(p):return json.loads(Path(p).read_text())
def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()

def run_r(command,script,args,data,log):
 cmd=[command,'--vanilla','--slave','-f',str(script),'--args']+list(map(str,args));t=time.perf_counter()
 with Path(log).open('w') as f:
  p=subprocess.run(cmd,input=data,text=True,stdout=f,stderr=subprocess.STDOUT)
 if p.returncode:raise RuntimeError(f'R failed with return code {p.returncode}; see {log}')
 return {'command':cmd,'elapsed_seconds':time.perf_counter()-t,'returncode':p.returncode}

def build_inputs(out):
 source=ROOT/'inputs/common_input_1573.csv';assert sha(source)==SOURCE_SHA
 rr=rows(source);by={r['insref_id']:r for r in rr};assert len(rr)==len(by)==1573;built={};checks=[]
 for key in KEYS:
  meta=load(ROOT/'inputs/supports'/key/'INPUT_METADATA.json');refs=meta['source_references'];ids=[r['ID'] for r in refs]
  assert ids==sorted(ids) and len(ids)==len(set(ids));X=[];y=[];cl=[]
  for ref in refs:
   r=by[ref['ID']];assert r['fiscal_year']=='2024';c=ref['prefecture_code'];cl.append(c);y.append(float(r['collection_ratio_fraction']))
   X.append([1.,float(ref['D'])]+[float(r[n]) for n in meta['columns'][2:6]]+[float(c==j) for j in range(2,meta['G']+1)])
  X=np.asarray(X,dtype='<f8');y=np.asarray(y,dtype='<f8');cl=np.asarray(cl,dtype='<i4');buf=io.StringIO(newline='');w=csv.writer(buf);w.writerow(['y']+meta['columns'][1:]+['pref'])
  for yy,xx,cc in zip(y,X,cl):w.writerow([repr(float(yy))]+[repr(float(v)) for v in xx[1:]]+[int(cc)])
  text=buf.getvalue();h={'X_sha256':bh(X.tobytes()),'y_sha256':bh(y.tobytes()),'cluster_int32_sha256':bh(cl.tobytes()),'stdin_csv_sha256':bh(text.encode())};test={n:h[n]==meta['memory_fingerprints'][n] for n in h};assert all(test.values()),(key,test)
  p=out/'inputs'/key;p.mkdir(parents=True);(p/'model_stdin.csv').write_text(text,newline='');(p/'EXPECTED_MEMORY_HASHES.txt').write_text('\n'.join(h[n] for n in ['X_sha256','y_sha256','cluster_int32_sha256'])+'\n');checks.append({'key':key,'checks':test,'hashes':h,'n':len(ids),'G':meta['G'],'p':X.shape[1]});built[key]={'X':X,'y':y,'cluster':cl,'ids':ids,'meta':meta,'hashes':h,'text':text}
 dump(out/'FIXED_INPUT_IDENTITY.json',{'source_sha256':SOURCE_SHA,'all_24_exact':True,'historical_raw_source_or_support_selection_builder_rerun':False,'checks':checks})
 overlap=sorted(set(built['P1_k1']['ids'])&set(built['P3_k1']['ids']));assert len(overlap)==15
 dump(out/'FIXED_SUPPORT_COUNTS_AND_OVERLAP.json',{'P1_P3_1SD_overlap_n':len(overlap),'P1_P3_overlap_IDs':overlap,'groups':[{'key':k,'n':b['meta']['n'],'G':b['meta']['G'],'original_D1':b['meta']['n_D1'],'original_D0':b['meta']['n_D0'],'display_first_n':b['meta']['n_D1'] if k.startswith('P1') else b['meta']['n_D0'],'display_second_n':b['meta']['n_D0'] if k.startswith('P1') else b['meta']['n_D1']} for k,b in built.items()]})
 return built

def reconstruct(out,built):
 sp=importlib.util.spec_from_file_location('independent_reconstruction',ROOT/'code/reconstruct.py');m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);records=[];composition=[]
 for key,b in built.items():
  X,y,cl,meta=b['X'],b['y'],b['cluster'],b['meta'];p=out/'points'/key;co=np.array([float(r['estimate']) for r in rows(p/'FRACTIONAL_COEFFICIENTS.csv')]);fd=rows(p/'FRACTIONAL_FIT_DIAGNOSTICS.csv')[0];X0=X.copy();X1=X.copy();X0[:,1]=0;X1[:,1]=1
  probs={n:np.array([m.expit(float(z)) for z in m.linpred(xx,co)]) for n,xx in [('observed',X),('D0',X0),('D1',X1)]};delta=100*math.fsum(float(z) for z in probs['D1']-probs['D0'])/len(y);score=np.array([math.fsum(float(x)*float(res) for x,res in zip(X[:,j],y-probs['observed'])) for j in range(X.shape[1])]);den=np.array([1+math.fsum(abs(float(x)) for x in X[:,j]) for j in range(X.shape[1])]);norm=float(np.max(abs(score)/den));rank=m.rk(X*np.sqrt(probs['observed']*(1-probs['observed']))[:,None])['rank']
  record={'key':key,'independent_delta_pp':delta,'original_R_delta_pp':float(fd['delta_pp']),'absolute_delta_difference_pp':abs(delta-float(fd['delta_pp'])),'independent_normalized_score':norm,'fixed_rank':m.rk(X)['rank'],'interior_rank':m.rk(X[(y>0)&(y<1)])['rank'],'final_weighted_rank':rank,'columns':X.shape[1],'fits_in_this_script':0,'pass':abs(delta-float(fd['delta_pp']))<=1e-6 and norm<=1e-8 and rank==X.shape[1]};assert record['pass'];records.append(record)
  for dg in [1,0]:
   mask=X[:,1]==dg;z={'contrast':key[:2],'support_SD':meta['support_caliper'],'group':'D'+str(dg),'insurer_count':int(mask.sum()),'prefecture_count':len(set(cl[mask]))}
   for j,n in enumerate(meta['columns'][2:6],2):
    v=X[mask,j]
    for label,value in [('mean',np.mean(v)),('SD',np.std(v,ddof=1)),('median',np.median(v)),('min',np.min(v)),('max',np.max(v))]:z[n+'_'+label]=float(value)
   composition.append(z)
 dump(out/'FRACTIONAL_RECONSTRUCTION_FRESH.json',records);dump(out/'SUPPORT_COMPOSITION_FRESH.json',composition)

def fit_points(out,command,built):
 reports=[]
 for key,b in built.items():
  p=out/'points'/key;p.mkdir(parents=True);h=b['hashes'];do_ols=key in ['P1_k0p5','P3_k0p5']
  rr=run_r(command,ROOT/'code/fit_points.R',[p,'public_run',h['X_sha256'],h['y_sha256'],h['cluster_int32_sha256'],str(do_ols).upper(),'TRUE'],b['text'],p/'R_CONSOLE.log')
  fs=['FRACTIONAL_COEFFICIENTS.csv','FRACTIONAL_FIT_DIAGNOSTICS.csv','FRACTIONAL_DIAGNOSTIC_ARRAYS.csv','FRACTIONAL_SCORE.csv','FRACTIONAL_PREDICTION_DIAGNOSTICS.csv','FRACTIONAL_DERIVED_DIAGNOSTICS.rds']+(['OLS_COEFFICIENTS.csv','OLS_FIT_DIAGNOSTICS.csv','OLS_PREDICTION_DIAGNOSTICS.csv','LEAVE_ONE_PREFECTURE_OUT.csv'] if do_ols else [])
  for f in fs:
   if not (p/f).is_file():raise RuntimeError(f'Missing fresh point output: {key}/{f}')
  fd=rows(p/'FRACTIONAL_FIT_DIAGNOSTICS.csv')[0]
  assert int(fd['formal_fit_calls'])==1 and int(fd['n'])==b['meta']['n']
  reports.append({'key':key,**rr,'fresh_GLM_calls':1,'fresh_OLS_calls':int(do_ols),'output_sha256':{f:sha(p/f) for f in fs}})
  print('FRESH POINT COMPLETE',key,flush=True)
 dump(out/'POINT_RUN_SUMMARY.json',reports)
 return reports

def summarize_boundaries(p,label):
 records=rows(p/'BOUNDARIES.csv')
 expected={f'{label}_j{j:03d}_{side}' for j in range(1,201) for side in ['lower','upper']}
 assert len(records)==400 and {r['id'] for r in records}==expected
 by={(int(r['j']),r['side']):r for r in records};widths=[]
 for j in range(1,201):
  pair=[by[j,s] for s in ['lower','upper']]
  held=any(r['diagnostic_hold']=='TRUE' for r in pair)
  raw=all(r['official_p_strictly_below_alpha']=='TRUE' for r in pair)
  status='HOLD' if held else 'PASS' if raw else 'NONPASS'
  widths.append({'contrast':label,'j':j,'delta_pp':j/100,'tost_status':status})
 with (p/'THRESHOLDS.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(widths[0]));w.writeheader();w.writerows(widths)
 return {'contrast':label,'fresh_boundaries':len(records),'widths':len(widths),
  'width_counts':{s:sum(r['tost_status']==s for r in widths) for s in ['PASS','NONPASS','HOLD']},
  'first_pass_grid_pp':next((r['delta_pp'] for r in widths if r['tost_status']=='PASS'),None)}

def boundaries(out,command,built):
 reports=[]
 for label in ['P1','P2','P3','P2_k0p5']:
  half=label=='P2_k0p5';key=label if half else label+'_k1';b=built[key];p=out/'tost'/label;p.mkdir(parents=True)
  script=ROOT/'code'/('run_boundaries_half.R' if half else 'run_boundaries_baseline.R')
  (p/'EXPECTED_MEMORY_HASHES.txt').write_text('\n'.join(b['hashes'][n] for n in ['X_sha256','y_sha256','cluster_int32_sha256'])+'\n')
  run_r(command,ROOT/'code/verify_stdin_roundtrip.R',[p],b['text'],p/'ROUNDTRIP.log')
  for name,field in [('X_LE_F64_C','X_sha256'),('Y_LE_F64','y_sha256'),('CLUSTER_LE_I32','cluster_int32_sha256')]:
   assert (p/f'R_STDIN_{name}.sha256').read_text().split()[0]==b['hashes'][field]
  execution_key=bh((sha(script)+b['hashes']['stdin_csv_sha256']).encode())
  args=[p,label,'400']+([execution_key] if half else [])
  print('FRESH FULL TOST START',label,now(),flush=True)
  info=run_r(command,script,args,b['text'],p/'R_CONSOLE.log')
  report={**summarize_boundaries(p,label),**info,'seed':20261015,'requested_B':99999}
  reports.append(report);dump(out/'TOST_RUN_SUMMARY.json',reports)
  print('FRESH FULL TOST COMPLETE',label,report['elapsed_seconds'],flush=True)
 return reports

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--out',type=Path,required=True,help='New directory for all fresh results')
 ap.add_argument('--r-command',required=True,help='Pinned R executable or launcher, not Rscript')
 ap.add_argument('--points-only',action='store_true',help='Partial smoke: six fractional fits and two half-SD OLS fits; no bootstrap')
 args=ap.parse_args();out=args.out.resolve()
 if out.exists():raise SystemExit('--out must be a new directory; existing outputs are never reused')
 if out==ROOT or ROOT in out.parents:raise SystemExit('--out must be outside the source package')
 command=shutil.which(args.r_command) or args.r_command
 command=str(Path(command).resolve())
 if not Path(command).is_file():raise SystemExit('R executable or launcher not found')
 out.mkdir(parents=True)
 receipt={'started_utc':now(),'mode':'points_only' if args.points_only else 'full_1600_boundaries','status':'RUNNING',
  'python_version':sys.version,'numpy_version':np.__version__,'support_selection_rerun':False}
 dump(out/'RUN_SUMMARY.json',receipt)
 try:
  expr="stopifnot(getRversion()=='4.6.1',packageVersion('stats')=='4.6.1',packageVersion('fwildclusterboot')=='0.14.3',packageVersion('dqrng')=='0.4.1');print(sessionInfo())"
  preflight=subprocess.run([command,'--vanilla','--slave','-e',expr],text=True,capture_output=True)
  (out/'R_ENVIRONMENT.txt').write_text(preflight.stdout+preflight.stderr)
  if preflight.returncode:raise RuntimeError('Pinned R/package runtime preflight failed; see R_ENVIRONMENT.txt')
  built=build_inputs(out);fit_points(out,command,built);reconstruct(out,built)
  tost=[] if args.points_only else boundaries(out,command,built)
  receipt.update({'completed_utc':now(),'status':'COMPLETE_PARTIAL_POINTS_ONLY' if args.points_only else 'COMPLETE_FULL_FIXED_SUPPORT',
   'fresh_fractional_models':6,'fresh_TOST_boundaries':sum(x['fresh_boundaries'] for x in tost),'fresh_TOST_widths':sum(x['widths'] for x in tost),
   'historical_result_comparison_performed':False})
  dump(out/'RUN_SUMMARY.json',receipt);print(json.dumps(receipt,indent=2))
 except BaseException as e:
  receipt.update({'status':'INCOMPLETE_OR_FAILED','stopped_utc':now(),'error':repr(e)})
  dump(out/'RUN_SUMMARY.json',receipt);raise
if __name__=='__main__':main()
