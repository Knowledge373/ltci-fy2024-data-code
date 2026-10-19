#!/usr/bin/env python3
"""Approved fixed-spec reconstruction, not recovered-original code replay."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
from pathlib import Path
import argparse,csv,datetime,hashlib,importlib.util,json,platform,sys,traceback
from decimal import Decimal,ROUND_HALF_EVEN
import numpy as np
from scipy import linalg,stats
import scipy
ROOT=Path(__file__).resolve().parent
LABEL='REIMPLEMENTED_FROM_DOCUMENTED_FIXED_SPECIFICATION'
CHECKS=[]
POLICY=['annual_4','annual_5','annual_6','annual_7','annual_9','annual_10','annual_12','provisional_none']
COV=['income_low3_share','age_75_84_share','age_85plus_share','ln_insured_end_fy_count']
ATOL,DF_ATOL,RTOL=1e-10,1e-8,1e-8

def now():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def rows(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def csvout(p,data):
    with Path(p).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def check(name,ok,detail=None):
    CHECKS.append({'check':name,'pass':bool(ok),'detail':detail})
    if not ok:raise RuntimeError('FAIL_STOP: '+name+' '+str(detail))
def compare(name,a,b,df=False):
    a,b=np.asarray(a,float),np.asarray(b,float)
    check(name+'_shape',a.shape==b.shape)
    tol=(DF_ATOL if df else ATOL)+RTOL*np.maximum(abs(a),abs(b));d=abs(a-b)
    check(name,np.isfinite(a).all() and np.isfinite(b).all() and (d<=tol).all(),{'max_abs_difference':float(d.max()),'max_tolerance_ratio':float((d/tol).max())})
def rank_check(name,X,expected):
    _,r,_=linalg.qr(X,mode='economic',pivoting=True);d=abs(np.diag(r));v={str(t):int((d>t*d.max()).sum()) for t in [1e-7,1e-9,1e-11]}
    check(name,all(n==expected for n in v.values()),v);return v

def load_module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def render(v,digits=3,signed=False):
    q=Decimal(str(float(v))).quantize(Decimal(1).scaleb(-digits),rounding=ROUND_HALF_EVEN)
    return ('+' if signed and q>=0 else '')+format(q,f'.{digits}f')

def table_expected():
    return {}

def add_display(display, expected, table, row, col, actual, raw=None):
    pass

def make_cell_design(annual,none,pref,C,reference):
    cells=sorted(set(zip(annual.tolist(),none.tolist())))
    nonref=[x for x in cells if x!=reference]
    names=['intercept']+[f'cell_{a}_{d}' for a,d in nonref]+COV+[f'prefecture_{g:02d}' for g in range(2,48)]
    X=np.column_stack([np.ones(len(pref))]+[((annual==a)&(none==d)).astype(float) for a,d in nonref]+[C[:,j] for j in range(4)]+[(pref==g).astype(float) for g in range(2,48)])
    return X,names,cells

def vector(names,positive,negative,reference):
    c=np.zeros(len(names))
    for cell,sign in [(positive,1),(negative,-1)]:
        if cell!=reference:c[names.index(f'cell_{cell[0]}_{cell[1]}')]=sign
    return c

def direct_deletion_fit(keep,T,W,y,label):
    # Every deletion is directly refitted. Exactly absent nuisance columns only.
    t=T[keep];active=np.any(t!=0,axis=0);t=t[:,active];w=W[keep];yy=y[keep]
    rank_check(label+'_nuisance_rank',t,t.shape[1])
    Q,_=linalg.qr(t,mode='economic');wr=w-Q@(Q.T@w);yr=yy-Q@(Q.T@yy)
    rank_check(label+'_target_rank',wr,12)
    q,r=linalg.qr(wr,mode='economic');b=linalg.solve_triangular(r,q.T@yr)
    # Independent direct full-span QR, with target columns first.
    z=np.column_stack([w,t]);q2,r2=linalg.qr(z,mode='economic');b2=linalg.solve_triangular(r2,q2.T@yy)[:12]
    compare(label+'_partial_vs_direct_full',b,b2)
    return b,int(t.shape[1])

def main(out):
    check('runtime_exact',platform.python_version()=='3.12.14' and np.__version__=='2.3.5' and scipy.__version__=='1.17.0')
    fixed={str(p.relative_to(ROOT)):sha(p) for folder in ['inputs','original_algebra','portable_algebra'] for p in (ROOT/folder).rglob('*') if p.is_file()}
    for p in [Path(__file__)]:fixed[str(p.relative_to(ROOT))]=sha(p)
    save(out/'PRE_EXECUTION_FREEZE.json',{'status':'FROZEN_BEFORE_RECONSTRUCTION','utc':now(),'label':LABEL,'hashes':fixed,'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},'all_direct_deletions':True,'high_hat_threshold_invented':False,'original_unrounded_results_recovered':False})
    rr=rows(ROOT/'inputs/common_input_1573.csv');spec=json.loads((ROOT/'inputs/baseline_model_specification.json').read_text())
    check('csv_pinned',sha(ROOT/'inputs/common_input_1573.csv')=='ea25bf07ed1ff42626ca73719ca252af7f563ca35cd19029139863cdbca468f2')
    ids=[r['insref_id'] for r in rr];annual=np.array([int(r['annual_installments']) for r in rr]);none=np.array([int(r['provisional_none']) for r in rr]);pref=np.array([int(r['prefecture_code']) for r in rr],dtype='<i4')
    C=np.array([[float(r[c]) for c in COV] for r in rr]);y=np.array([float(r['collection_ratio_fraction']) for r in rr],dtype='<f8');names61=spec['columns']
    X61=np.asarray([[1.0]+[float(int(r['annual_installments'])==n) for n in range(3,13) if n!=8]+[float(r['provisional_none'])]+[float(r[c]) for c in COV]+[float(int(r['prefecture_code'])==n) for n in range(2,48)] for r in rr],dtype='<f8')
    check('baseline_X_hash',hashlib.sha256(X61.tobytes()).hexdigest()=='5cf2c6bfb14e920e6889ed9c690a3eb21c359eb130ddfe9bf00b3e4c171144d7')
    check('outcome_y_hash',hashlib.sha256(y.tobytes()).hexdigest()=='b078ce93f589e269062f8784526e70fb2a26a20ac354382d577b944b5bc98719')
    check('cluster_hash',hashlib.sha256(pref.tobytes()).hexdigest()=='3ae9f549001855930cff396274e7adfcda7e9a331b84ac2f709b0468d694893c')
    check('row_ids',len(ids)==1573 and ids==sorted(ids) and len(set(ids))==1573)
    rank_check('baseline_full61',X61,61);qb,rb=linalg.qr(X61,mode='economic');beta61=linalg.solve_triangular(rb,qb.T@y);fit61=qb@(qb.T@y)
    expected=table_expected();display=[]
    baseline_diag={'model':'national_additive61','n':len(y),'p':61,'minimum_predicted_percent':float(fit61.min()*100),'maximum_predicted_percent':float(fit61.max()*100),'below_zero':int((fit61<0).sum()),'above_one':int((fit61>1).sum()),'rmse_pp':float(np.sqrt(np.mean((y-fit61)**2))*100)}
    for col,key in [(2,'minimum_predicted_percent'),(3,'maximum_predicted_percent')]:add_display(display,expected,'S9',1,col,render(baseline_diag[key]),baseline_diag[key])
    add_display(display,expected,'S9',1,4,str(baseline_diag['above_one']),baseline_diag['above_one'])
    save(out/'NATIONAL_PREDICTION_DIAGNOSTICS.json',baseline_diag)
    # Run unchanged CR2 algebra with documented fixed16 model.
    sys.path.insert(0,str(ROOT/'portable_algebra'))
    primary=load_module('existing_cr2_algebra',ROOT/'portable_algebra/run_two_contrast_cr2.py')
    independent=load_module('existing_independent_algebra',ROOT/'portable_algebra/run_independent_projector_qa.py')
    (out/'independent_cell16').mkdir();independent.OUT=out/'independent_cell16'
    X,names,cells=make_cell_design(annual,none,pref,C,(8,1));check('sixteen_cells',len(cells)==16)
    singleton=[cell for cell in cells if int(((annual==cell[0])&(none==cell[1])).sum())==1]
    check('singleton_cells',singleton==[(3,1),(7,0),(11,0)],singleton);rank_check('cell_full66',X,66)
    T=np.column_stack([(pref==g).astype(float) for g in range(1,48)]+[((annual==a)&(none==d)).astype(float) for a,d in singleton]);rank_check('cell_nuisance50',T,50)
    target_names=[f'cell_{a}_{d}' for a,d in cells if (a,d)!=(8,1) and (a,d) not in singleton]+COV;idx=[names.index(t) for t in target_names]
    Qt,_=linalg.qr(T,mode='economic');P=X[:,idx]-Qt@(Qt.T@X[:,idx]);yp=y-Qt@(Qt.T@y)
    for g in range(48):rank_check('cell_partial16_omit_pref_'+str(g),P if g==0 else P[pref!=g],16)
    smask=np.array([(a,d) in singleton for a,d in zip(annual,none)])
    compare('cell_singleton_partial_X',P[smask],np.zeros_like(P[smask]));compare('cell_singleton_partial_y',yp[smask],np.zeros_like(yp[smask]))
    full=primary.calc_model(X,y,pref,names,'full66');part=primary.calc_model(P,yp,pref,target_names,'partial16')
    check('cell_full_suppressed50',sum(r['suppressed_1e12'] for r in full['eig'])==50)
    check('cell_partial_suppressed0',sum(r['suppressed_1e12'] for r in part['eig'])==0)
    check('cell_no_material_negative_eigen',all(r['material_negative_count']==0 for r in full['eig']+part['eig']))
    compare('cell_full_partial_target_beta',full['beta'][idx],part['beta']);compare('cell_full_partial_residuals',full['e'],part['e']);compare('cell_full_partial_covariance',full['V'][np.ix_(idx,idx)],part['V'])
    ifit=independent.fit_representation(X,y,pref,'cell_original66')
    altX,alt_names,_=make_cell_design(annual,none,pref,C,(6,0));alt=primary.calc_model(altX,y,pref,alt_names,'alternate_ref_full66')
    compare('cell_alternative_reference_residuals',full['e'],alt['e'])
    result=[];contrasts=[('C1',(6,1),(6,0)),('C2',(10,1),(10,0)),('C3',(9,1),(8,1)),('C4',(12,0),(6,0))]
    for j,(label,pos,neg) in enumerate(contrasts,1):
        cv=vector(names,pos,neg,(8,1));a=primary.contrast(full,cv);b=primary.contrast(part,cv[idx]);d,_=independent.infer_from_projector(ifit,cv,y);other=primary.contrast(alt,vector(alt_names,pos,neg,(6,0)))
        for field,df in [('estimate_fraction',False),('variance_fraction_squared',False),('SE_fraction',False),('df',True),('CI_lower_fraction',False),('CI_upper_fraction',False)]:
            compare(label+'_full_partial_'+field,a[field],b[field],df);compare(label+'_reference_invariant_'+field,a[field],other[field],df);indkey='variance_fraction2' if field=='variance_fraction_squared' else field;compare(label+'_independent_'+field,a[field],d[indkey],df)
        compare(label+'_GH_vs_direct_P',a['P'],a['Pdirect']);compare(label+'_full_partial_P',a['P'],b['P'])
        posmask=(annual==pos[0])&(none==pos[1]);negmask=(annual==neg[0])&(none==neg[1]);both=sorted(set(pref[posmask])&set(pref[negmask]));within=np.isin(pref,both)
        record={'contrast':label,'positive_cell':f'{pos[0]}_'+str(pos[1]),'negative_cell':f'{neg[0]}_'+str(neg[1]),'estimate_fraction':a['estimate_fraction'],'estimate_pp':100*a['estimate_fraction'],'variance_fraction_squared':a['variance_fraction_squared'],'SE_fraction':a['SE_fraction'],'SE_pp':100*a['SE_fraction'],'df':a['df'],'t975':a['t975'],'CI_lower_pp':100*a['CI_lower_fraction'],'CI_upper_pp':100*a['CI_upper_fraction'],'positive_n':int(posmask.sum()),'negative_n':int(negmask.sum()),'positive_prefectures':len(set(pref[posmask])),'negative_prefectures':len(set(pref[negmask])),'overlap_prefectures':len(both),'positive_overlap_n':int((posmask&within).sum()),'negative_overlap_n':int((negmask&within).sum()),'provenance':LABEL}
        result.append(record)
        for col,key in [(1,'estimate_pp'),(2,'CI_lower_pp'),(3,'CI_upper_pp'),(4,'df')]:add_display(display,expected,'S4',j,col,render(record[key],signed=col in [1,2,3]),record[key])
        add_display(display,expected,'S4',j,5,str(record['overlap_prefectures']),record['overlap_prefectures'])
    csvout(out/'CELL16_FOUR_CONTRASTS.csv',result);save(out/'CELL16_FOUR_CONTRASTS.json',result);csvout(out/'CELL16_EIGEN_DIAGNOSTICS.csv',full['eig']+part['eig'])
    fitted=y-full['e'];cell_diag={'model':'documented_cell16_66columns','n':len(y),'p':66,'minimum_predicted_percent':float(fitted.min()*100),'maximum_predicted_percent':float(fitted.max()*100),'below_zero':int((fitted<0).sum()),'above_one':int((fitted>1).sum()),'rmse_pp':float(np.sqrt(np.mean(full['e']**2))*100)}
    for col,key in [(2,'minimum_predicted_percent'),(3,'maximum_predicted_percent')]:add_display(display,expected,'S9',8,col,render(cell_diag[key]),cell_diag[key])
    add_display(display,expected,'S9',8,4,str(cell_diag['above_one']),cell_diag['above_one']);save(out/'CELL16_PREDICTION_DIAGNOSTICS.json',cell_diag)
    check('cell_rmse_display',render(cell_diag['rmse_pp'])=='2.296',cell_diag['rmse_pp'])
    print('Fixed cell-model calculations complete; starting all-direct deletion checks.',flush=True)
    save(out/'CHECKS_INTERIM_CELL16.json',CHECKS);save(out/'DISPLAY_MATCHES_INTERIM_CELL16.json',display)
    # Existing additive61 partial leverage and deletion diagnostic definitions.
    influences=[];levs=[];pref_shares=[];cell_shares=[]
    for term in POLICY:
        k=names61.index(term);qo,_=linalg.qr(np.delete(X61,k,axis=1),mode='economic');r=X61[:,k]-qo@(qo.T@X61[:,k]);v=r*r/(r@r)
        compare('leverage_'+term+'_sum',sum(v),1)
        ps={g:float(v[pref==g].sum()) for g in range(1,48)};cs={cell:float(v[(annual==cell[0])&(none==cell[1])].sum()) for cell in cells}
        compare('leverage_'+term+'_pref_sum',sum(ps.values()),1);compare('leverage_'+term+'_cell_sum',sum(cs.values()),1)
        mxp=max(ps.values());mxc=max(cs.values());pref_ties=[g for g in ps if ps[g]==mxp];cell_ties=[c for c in cs if cs[c]==mxc]
        influences.append({'term':term,'baseline_coefficient_pp':float(beta61[k]*100),'maximum_prefecture_share_percent':mxp*100,'maximum_prefecture_codes':pref_ties,'maximum_prefecture_names':[spec['prefecture_codes'][f'{g:02d}'] for g in pref_ties],'maximum_cell_share_percent':mxc*100,'maximum_cells':cell_ties})
        levs.extend({'insref_id':i,'term':term,'partial_leverage':float(x)} for i,x in zip(ids,v));pref_shares.extend({'term':term,'prefecture':g,'partial_leverage':value} for g,value in ps.items());cell_shares.extend({'term':term,'annual':c[0],'provisional_none':c[1],'partial_leverage':value} for c,value in cs.items())
    widx=[names61.index(x) for x in POLICY+COV];W=X61[:,widx];T0=np.column_stack([(pref==g).astype(float) for g in range(1,48)]+[(annual==a).astype(float) for a in [3,11]])
    base,_=direct_deletion_fit(np.ones(len(y),bool),T0,W,y,'S3_baseline');compare('S3_baseline_original',base,beta61[widx]);deletions=[];pref_b=[];ins_b=[]
    for g in range(1,48):
        b,rank=direct_deletion_fit(pref!=g,T0,W,y,'S3_delete_pref_'+str(g));pref_b.append(b)
        deletions.append({'deleted_type':'prefecture','deleted_id':f'{g:02d}','retained_n':int((pref!=g).sum()),'nuisance_rank':rank,**{term:float(b[k]*100) for k,term in enumerate(POLICY)}})
    for i,identity in enumerate(ids):
        keep=np.arange(len(y))!=i;b,rank=direct_deletion_fit(keep,T0,W,y,'S3_delete_insurer_'+identity);ins_b.append(b)
        deletions.append({'deleted_type':'insurer','deleted_id':identity,'retained_n':len(y)-1,'nuisance_rank':rank,**{term:float(b[k]*100) for k,term in enumerate(POLICY)}})
        if (i+1)%400==0:print('Direct insurer deletions checked:',i+1,flush=True)
    pref_b=np.array(pref_b);ins_b=np.array(ins_b)
    for j,record in enumerate(influences):
        mn,mx=float(pref_b[:,j].min()*100),float(pref_b[:,j].max()*100);delta=abs(ins_b[:,j]-base[j])*100;maxdelta=float(delta.max());ties=[ids[k] for k in np.flatnonzero(delta==delta.max())]
        record.update({'prefecture_deletion_min_pp':mn,'prefecture_deletion_max_pp':mx,'maximum_insurer_absolute_change_pp':maxdelta,'maximum_insurer_ids':ties,'all_deletions_directly_refitted':True,'provenance':LABEL})
        check('S3_'+record['term']+'_unique_pref_cell_max',len(record['maximum_prefecture_codes'])==len(record['maximum_cells'])==1)
        cell=record['maximum_cells'][0];cell_label=f'{cell[0]}回・'+('なし' if cell[1] else 'あり')
        add_display(display,expected,'S3',j+1,1,render(record['maximum_prefecture_share_percent'],2)+'（'+record['maximum_prefecture_names'][0]+'）',record['maximum_prefecture_share_percent'])
        add_display(display,expected,'S3',j+1,2,cell_label+'／'+render(record['maximum_cell_share_percent'],2),record['maximum_cell_share_percent'])
        add_display(display,expected,'S3',j+1,3,render(mn,signed=True)+'～'+render(mx,signed=True),{'minimum_pp':mn,'maximum_pp':mx})
        add_display(display,expected,'S3',j+1,4,render(maxdelta),maxdelta)
    save(out/'S3_INFLUENCE_SUMMARY.json',influences);csvout(out/'S3_ALL_DELETION_COEFFICIENTS_PP.csv',deletions);csvout(out/'S3_INSURER_PARTIAL_LEVERAGE.csv',levs);csvout(out/'S3_PREFECTURE_PARTIAL_LEVERAGE.csv',pref_shares);csvout(out/'S3_CELL_PARTIAL_LEVERAGE.csv',cell_shares)
    save(out/'DISPLAY_MATCHES.json',display);save(out/'CHECKS.json',CHECKS)
    after={name:sha(ROOT/name) for name in fixed};check('protected_sources_unchanged',after==fixed)
    check('independent_cell16_all_checks_pass',all(x['pass'] for x in independent.CHECKS))
    save(out/'CHECKS.json',CHECKS)
    receipt={'status':'PASS_FIXED_SPEC_RECONSTRUCTION','label':LABEL,'completed_utc':now(),'checks':len(CHECKS),'independent_cell16_checks':len(independent.CHECKS),'s3_direct_prefecture_deletions':47,'s3_direct_insurer_deletions':1573,'original_code_recovered':False,'historical_unrounded_equivalence_claimed':False,'manuscript_comparison_performed':False,'fresh_precise_values_verified_by_independent_routes':True,'model_changed':False,'classification_changed':False}
    save(out/'RECONSTRUCTION_RECEIPT.json',receipt);print(json.dumps(receipt,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output-dir',type=Path);a=p.parse_args();out=(a.output_dir or ROOT/'runs'/datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')).resolve()
    if out.exists():raise SystemExit('Refusing existing output directory.')
    out.mkdir(parents=True)
    try:main(out)
    except Exception as exc:
        save(out/'CHECKS_BEFORE_FAILURE.json',CHECKS);save(out/'FAILURE.json',{'utc':now(),'status':'FAIL_STOP','error':repr(exc),'traceback':traceback.format_exc(),'label':LABEL});raise
