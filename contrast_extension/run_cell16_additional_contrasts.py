#!/usr/bin/env python3
"""Newly rebuilt thin runner; old run001 was not recovered. Separate new run only."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
from pathlib import Path
import argparse, csv, datetime, hashlib, importlib.util, json, platform, traceback, warnings
import numpy as np
import scipy
from scipy import linalg

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parent
CHECKS = []
FIELDS = ['estimate_fraction','variance_fraction_squared','SE_fraction','df','t975','CI_lower_fraction','CI_upper_fraction']
COUNT_FIELDS = ['positive_n','negative_n','positive_prefectures','negative_prefectures','overlap_prefectures','positive_overlap_n','negative_overlap_n']
PROVENANCE = 'NEW_REBUILT_EXECUTION_OF_FIXED_EXPLORATORY_CELL16_CONTRASTS_PRIOR_RUN_NOT_RECOVERED'

def now(): return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p, data): Path(p).write_text(json.dumps(data,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def rows(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def csvout(p, data):
    with Path(p).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def check(name, ok, detail=None):
    CHECKS.append({'check':name,'pass':bool(ok),'detail':detail})
    if not ok: raise RuntimeError('HOLD: '+name+' '+str(detail))
def compare(name, a, b, df=False):
    a,b=np.asarray(a,float),np.asarray(b,float)
    check(name+'_shape',a.shape==b.shape)
    tol=(1e-8 if df else 1e-10)+1e-8*np.maximum(abs(a),abs(b));d=abs(a-b)
    check(name,a.size>0 and np.isfinite(a).all() and np.isfinite(b).all() and (d<=tol).all(),{'max_absolute_difference':float(d.max()),'max_tolerance_ratio':float((d/tol).max())})
def load(name, path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
def tree(path):
    return {str(p.relative_to(path)):sha(p) for p in sorted(Path(path).rglob('*')) if p.is_file()}
def protected(freeze):
    if freeze:
        for root,expected in freeze['protected_trees'].items(): check('protected_tree_'+root,tree(Path(root))==expected)
        for path,expected in freeze['protected_files'].items(): check('protected_file_'+path,Path(path).is_file() and sha(path)==expected)
def counts(annual,none,pref,pos,neg):
    a=(annual==pos[0])&(none==pos[1]);b=(annual==neg[0])&(none==neg[1])
    pa,pb=set(pref[a]),set(pref[b]);both=pa&pb;within=np.isin(pref,sorted(both))
    return dict(zip(COUNT_FIELDS,[int(a.sum()),int(b.sum()),len(pa),len(pb),len(both),int((a&within).sum()),int((b&within).sum())]))
def numeric(a):
    return {**{k:a[k] for k in FIELDS},'estimate_pp':100*a['estimate_fraction'],'SE_pp':100*a['SE_fraction'],'CI_lower_pp':100*a['CI_lower_fraction'],'CI_upper_pp':100*a['CI_upper_fraction']}

def run(out, freeze, freeze_path, prior_path):
    spec=json.loads((HERE/'SPECIFICATION.json').read_text())
    versions={'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'executable':sys.executable}
    check('runtime_supported',platform.python_version_tuple()[:2]==('3','12') and np.__version__=='2.3.5' and scipy.__version__=='1.17.0',versions)
    if versions['python']!='3.12.14': warnings.warn('Python patch differs from 3.12.14. Numerical checks remain mandatory.')
    before=tree(PACKAGE)
    protected(freeze)
    if freeze:
        check('preflight_source_hashes',before==freeze['package_file_hashes'])
        check('preflight_runtime',versions==freeze['runtime'])
        check('preflight_new_output',str(out)==freeze['new_empty_output_directory'])
    for path,h in spec['legacy_file_hashes'].items(): check('legacy_pinned_'+path,sha(PACKAGE/path)==h)
    for path,h in spec['reference_hashes'].items(): check('reference_pinned_'+path,sha(HERE/'reference'/path)==h)
    save(out/'PRE_EXECUTION_FREEZE.json',{'status':'FROZEN_BEFORE_NEW_REBUILT_MODEL_FIT','utc':now(),'specification_id':spec['specification_id'],'package_file_hashes':before,'runtime':versions,'external_freeze_sha256':sha(freeze_path) if freeze_path else None,'new_output_was_empty':True,'prior_run001_recovered':False,'prior_results_known':True,'new_preregistration_claimed':False,'baseline_public_commit':spec['baseline_public_commit'],'prior_results_transcription_sha256':sha(prior_path) if prior_path else None})
    save(out/'RUNTIME_ENVIRONMENT.json',versions)
    root=PACKAGE/'fixed_model_reconstruction'
    recon=load('unchanged_reconstruction',root/'run_reconstruction.py')
    primary=load('unchanged_cr2',root/'portable_algebra/run_two_contrast_cr2.py')
    independent=load('unchanged_projector',root/'portable_algebra/run_independent_projector_qa.py')
    # First independent function invocation is strictly AFTER this redirection.
    (out/'independent_cell16').mkdir()
    independent.OUT=out/'independent_cell16'
    check('independent_OUT_before_first_call',independent.OUT.resolve()==(out/'independent_cell16').resolve())
    rr=rows(PACKAGE/spec['input']);ids=[r['insref_id'] for r in rr]
    check('1573_unique_sorted_rows',len(ids)==1573 and len(set(ids))==1573 and ids==sorted(ids))
    annual=np.array([int(r['annual_installments']) for r in rr]);none=np.array([int(r['provisional_none']) for r in rr])
    pref=np.array([int(r['prefecture_code']) for r in rr],dtype='<i4')
    y=np.array([float(r['collection_ratio_fraction']) for r in rr],dtype='<f8')
    C=np.array([[float(r[k]) for k in recon.COV] for r in rr])
    check('y_memory_hash',hashlib.sha256(y.tobytes()).hexdigest()==spec['y_hash'])
    check('pref_memory_hash',hashlib.sha256(pref.tobytes()).hexdigest()==spec['pref_hash'])
    check('all_47_prefectures',sorted(set(pref.tolist()))==list(range(1,48)))
    X,names,cells=recon.make_cell_design(annual,none,pref,C,(8,1))
    check('exact_66_columns',names==spec['columns'] and X.shape==(1573,66))
    check('exact_16_cells',[list(c) for c in cells]==spec['observed_cells'])
    singleton=[c for c in cells if int(((annual==c[0])&(none==c[1])).sum())==1]
    check('singleton_cells_retained',singleton==[(3,1),(7,0),(11,0)])
    recon.rank_check('full66',X,66)
    T=np.column_stack([(pref==g).astype(float) for g in range(1,48)]+[((annual==a)&(none==d)).astype(float) for a,d in singleton])
    recon.rank_check('nuisance50',T,50)
    target_names=[f'cell_{a}_{d}' for a,d in cells if (a,d)!=(8,1) and (a,d) not in singleton]+recon.COV
    idx=[names.index(n) for n in target_names];Qt,_=linalg.qr(T,mode='economic')
    XP=X[:,idx]-Qt@(Qt.T@X[:,idx]);yp=y-Qt@(Qt.T@y)
    recon.rank_check('partial16',XP,16)
    smask=np.array([(a,d) in singleton for a,d in zip(annual,none)])
    compare('singleton_partial_X_zero',XP[smask],np.zeros_like(XP[smask]));compare('singleton_partial_y_zero',yp[smask],np.zeros_like(yp[smask]))
    full=primary.calc_model(X,y,pref,names,'full66');part=primary.calc_model(XP,yp,pref,target_names,'partial16')
    check('full_suppressed50',sum(r['suppressed_1e12'] for r in full['eig'])==50)
    check('partial_suppressed0',sum(r['suppressed_1e12'] for r in part['eig'])==0)
    check('no_material_negative_eigen',all(r['material_negative_count']==0 for r in full['eig']+part['eig']))
    compare('full_partial_target_beta',full['beta'][idx],part['beta'])
    compare('full_partial_residuals',full['e'],part['e']);compare('full_partial_covariance',full['V'][np.ix_(idx,idx)],part['V'])
    ifit=independent.fit_representation(X,y,pref,'cell_original66')
    historical=json.loads(prior_path.read_text())['rows'] if prior_path else None
    if historical is not None: check('prior_reference_ids',[r['contrast'] for r in historical]==['B9_cell16','B10_cell16'])
    result=[];routes=[];historical_checks=[];matrix_checks=[]
    for j,definition in enumerate(spec['new_contrasts']):
        old=historical[j] if historical is not None else None
        label=definition['id'];pos=tuple(definition['positive']);neg=tuple(definition['negative'])
        c=recon.vector(names,pos,neg,(8,1))
        check(label+'_indices_signs',{str(i):float(c[i]) for i in np.flatnonzero(c)}==definition['weights_by_zero_based_index'] and np.count_nonzero(c)==2 and c.sum()==0)
        check(label+'_named_signs',{names[i]:float(c[i]) for i in np.flatnonzero(c)}=={f'cell_{pos[0]}_{pos[1]}':1,f'cell_{neg[0]}_{neg[1]}':-1})
        a=primary.contrast(full,c);b=primary.contrast(part,c[idx]);d,scores=independent.infer_from_projector(ifit,c,y)
        d={**d,'variance_fraction_squared':d['variance_fraction2'],'t975':d['critical_t']}
        check(label+'_finite_positive',all(np.isfinite(a[k]) for k in FIELDS) and min(a['variance_fraction_squared'],a['SE_fraction'],a['df'])>0)
        for field in FIELDS:
            compare(label+'_full_partial_'+field,a[field],b[field],field=='df');compare(label+'_independent_'+field,a[field],d[field],field=='df')
            if old is not None:
                compare(label+'_prior_run_'+field,a[field],old[field],field=='df')
                historical_checks.append({'contrast':label,'field':field,'new_value':a[field],'prior_conversation_value':old[field],'absolute_difference':abs(a[field]-old[field]),'pass':True})
        compare(label+'_GH_direct_P',a['P'],a['Pdirect']);compare(label+'_full_partial_P',a['P'],b['P'])
        hf=np.stack([c@h for h in full['H']],axis=1);hp=np.stack([c[idx]@h for h in part['H']],axis=1)
        check(label+'_H_P_orientations',hf.shape==(66,47) and hp.shape==(16,47) and a['P'].shape==(47,47))
        matrix_checks.append({'contrast':label,'H_full_shape':list(hf.shape),'H_partial_shape':list(hp.shape),'P_shape':list(a['P'].shape),'qTy_score_identity':'PASS in unchanged independent infer_from_projector'})
        placement=counts(annual,none,pref,pos,neg)
        if old is not None:
            for field in COUNT_FIELDS:check(label+'_prior_placement_'+field,placement[field]==old[field])
        for route,value in [('full66',a),('partial16',b),('independent_projector',d)]: routes.append({'contrast':label,'route':route,**numeric(value)})
        r={'contrast':label,'label_ja':f'{pos[0]}回なし−{neg[0]}回あり','positive_cell':f'{pos[0]}_{pos[1]}','negative_cell':f'{neg[0]}_{neg[1]}','contrast_weights_by_name':json.dumps({names[i]:float(c[i]) for i in np.flatnonzero(c)},sort_keys=True),**numeric(a),**placement,'provenance':PROVENANCE,'status':'NEW_COMPUTATION_PENDING_SEPARATE_QA'}
        result.append(r);print(json.dumps(r,ensure_ascii=False),flush=True)
    check('both_new_rows_in_order',[r['contrast'] for r in result]==['B9_cell16','B10_cell16'])
    refs=rows(HERE/'reference/CELL16_FOUR_CONFIRMED.csv');check('old4_ids',[r['contrast'] for r in refs]==['C1','C2','C3','C4'])
    old_checks=[]
    definitions=[('C1',(6,1),(6,0)),('C2',(10,1),(10,0)),('C3',(9,1),(8,1)),('C4',(12,0),(6,0))]
    for (label,pos,neg),ref in zip(definitions,refs):
        a=primary.contrast(full,recon.vector(names,pos,neg,(8,1)));placement=counts(annual,none,pref,pos,neg)
        check(label+'_old_labels',ref['contrast']==label and ref['positive_cell']==f'{pos[0]}_{pos[1]}' and ref['negative_cell']==f'{neg[0]}_{neg[1]}')
        for field in FIELDS:
            value=float(ref[field]) if field in ref else float(ref[field.replace('_fraction','_pp')])/100
            compare(label+'_old_'+field,a[field],value,field=='df')
            old_checks.append({'contrast':label,'field':field,'new_value':a[field],'legacy_value':value,'absolute_difference':abs(a[field]-value),'pass':True})
        for field in COUNT_FIELDS:
            check(label+'_old_'+field,placement[field]==int(ref[field]));old_checks.append({'contrast':label,'field':field,'new_value':placement[field],'legacy_value':int(ref[field]),'absolute_difference':0,'pass':True})
        for field in ['estimate_pp','SE_pp','CI_lower_pp','CI_upper_pp','df']:
            check(label+'_old_display3_'+field,recon.render(numeric(a)[field])==recon.render(float(ref[field])))
    additive={r['contrast']:r for r in rows(HERE/'reference/ADDITIVE_TWO_CONFIRMED.csv')}
    check('additive_reference_ids',set(additive)=={'9no_minus_12yes','10no_minus_12yes'})
    comparison=[]
    for r,old_id in zip(result,['9no_minus_12yes','10no_minus_12yes']):
        old=additive[old_id]
        for model,numbers,source in [('additive61_legacy_saved',{k:float(old[k]) for k in FIELDS},'reference/ADDITIVE_TWO_CONFIRMED.csv'),('cell16_full66_rebuilt',r,'CELL16_ADDITIONAL_TWO_CONTRASTS.csv')]:
            nr=numeric(numbers)
            comparison.append({'contrast':r['contrast'],'label_ja':r['label_ja'],'model':model,**nr,**{k+'_display3':recon.render(nr[k],signed=k in ['estimate_pp','CI_lower_pp','CI_upper_pp']) for k in ['estimate_pp','SE_pp','df','CI_lower_pp','CI_upper_pp']},'source':source,'additive_reestimated':False})
    check('both_models_four_rows',len(comparison)==4)
    csvout(out/'CELL16_ADDITIONAL_TWO_CONTRASTS.csv',result);save(out/'CELL16_ADDITIONAL_TWO_CONTRASTS.json',result)
    csvout(out/'NUMERICAL_ROUTE_COMPARISON.csv',routes)
    if historical_checks:csvout(out/'PRIOR_RUN_NUMERICAL_COMPARISON.csv',historical_checks)
    csvout(out/'OLD4_REGRESSION_CHECK.csv',old_checks);csvout(out/'ADDITIVE_CELL16_COMPARISON.csv',comparison)
    csvout(out/'CELL16_EIGEN_DIAGNOSTICS.csv',full['eig']+part['eig']);save(out/'MATRIX_ORIENTATION_AND_SCORE_CHECKS.json',matrix_checks)
    csvout(out/'PREFECTURE_THREE_CELL_COUNTS.csv',[{'prefecture_code':f'{g:02d}',**{f'cell_{a}_{d}':int(((pref==g)&(annual==a)&(none==d)).sum()) for a,d in [(9,1),(10,1),(12,0)]}} for g in range(1,48)])
    CHECKS.extend(recon.CHECKS);check('all_independent_checks_pass',all(r['pass'] for r in independent.CHECKS))
    check('source_package_unchanged',tree(PACKAGE)==before);protected(freeze)
    save(out/'CHECKS.json',CHECKS)
    save(out/'EXECUTION_RECEIPT.json',{'status':'PASS_NEW_REBUILT_EXECUTION_PENDING_SEPARATE_INDEPENDENT_QA','completed_utc':now(),'specification_id':spec['specification_id'],'new_contrast_ids':['B9_cell16','B10_cell16'],'checks':len(CHECKS),'independent_projector_checks':len(independent.CHECKS),'prior_run001_recovered':False,'new_execution_after_environment_loss':True,'prior_results_known_before_new_run':True,'old4_regression_pass':True,'prior_run_numerical_agreement':True if historical is not None else None,'full_partial_projector_agree':True,'protected_hashes_unchanged':True,'manuscript_comparison_performed':False,'external_independent_QA_performed':False,'legacy_main_called':False,'additive_reestimated':False,'full_scope_or_bootstrap_or_TOST_run':False,'classification_month_rows_weights_supports_changed':False,'official_R_package_run':False,'publication_authorized':False,'provenance':PROVENANCE})

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output-dir',type=Path,required=True);ap.add_argument('--freeze',type=Path);ap.add_argument('--prior-results',type=Path,help='Optional local conversation-transcribed prior results for historical comparison only; not an analysis input.')
    args=ap.parse_args();out=args.output_dir.resolve();freeze=json.loads(args.freeze.read_text()) if args.freeze else None
    if out==PACKAGE or PACKAGE in out.parents: raise SystemExit('HOLD: output must be outside source package')
    if out.exists() and (not out.is_dir() or any(out.iterdir())): raise SystemExit('HOLD: output must be new or empty')
    if freeze:
        for root in freeze['protected_trees']:
            p=Path(root).resolve()
            if out==p or p in out.parents: raise SystemExit('HOLD: output overlaps protected tree')
    out.mkdir(parents=True,exist_ok=True)
    try: run(out,freeze,args.freeze,args.prior_results)
    except Exception as exc:
        save(out/'CHECKS.json',CHECKS);save(out/'FAILURE.json',{'status':'TECHNICAL_HOLD_NEW_REBUILD','utc':now(),'error':str(exc),'traceback':traceback.format_exc(),'both_contrasts_remain_in_scope':True})
        if not (out/'CELL16_ADDITIONAL_TWO_CONTRASTS.json').exists():
            failed=[{'contrast':c,'status':'NOT_FINALIZED_TECHNICAL_HOLD','reason':str(exc)} for c in ['B9_cell16','B10_cell16']]
            save(out/'CELL16_ADDITIONAL_TWO_CONTRASTS.json',failed);csvout(out/'CELL16_ADDITIONAL_TWO_CONTRASTS.csv',failed)
        raise
    finally: save(out/'OUTPUT_HASHES.json',{str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='OUTPUT_HASHES.json'})
