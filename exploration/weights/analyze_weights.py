#!/usr/bin/env python3
"""Fixed four-weight point-estimate analysis of the nationwide dataset."""
from __future__ import annotations
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import argparse,csv,hashlib,json,platform,sys
from pathlib import Path
from datetime import datetime,timezone
from decimal import Decimal,localcontext
import numpy as np
import scipy
from scipy import linalg

HERE=Path(os.environ['EXPLORATION_OUTPUT_DIR'])
BASE=Path(__file__).resolve().parents[1]/'inputs'/'baseline'
COV=['income_low3_share','age_75_84_share','age_85plus_share','ln_insured_end_fy_count']
FIXED=['insref_id','fiscal_year','prefecture_code','annual_installments','provisional_none','collection_ratio_fraction']+COV
SCHEMES=['equal_insurer','ordinary_current_assessment','all_first_category_insured','proxy_N_times_A_over_T']


def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open(encoding='utf-8-sig',newline='') as f:
        rr=csv.DictReader(f); return rr.fieldnames,list(rr)
def save(name,rows):
    assert rows,name
    with (HERE/name).open('w',encoding='utf-8-sig',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
def dump(name,obj):
    (HERE/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,allow_nan=False)+'\n',encoding='utf-8')
def scalar(v): return float(v)
def wmean(a,w): return float(np.dot(a,w)/w.sum())
def wcov(a,b,w): return wmean((a-wmean(a,w))*(b-wmean(b,w)),w)
def wcorr(a,b,w):
    den=np.sqrt(wcov(a,a,w)*wcov(b,b,w))
    return wcov(a,b,w)/den if den else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--expected-sha256',required=True)
    ap.add_argument('--count-column',required=True)
    ap.add_argument('--ordinary-assessment-column',required=True)
    ap.add_argument('--total-assessment-column',required=True)
    ap.add_argument('--ordinary-collection-column')
    args=ap.parse_args()
    started=datetime.now(timezone.utc).isoformat()
    assert sha(args.input)==args.expected_sha256,'Released input hash mismatch'
    fields,rows=read(args.input);bf,baseline=read(BASE/'analysis_rows_1573.csv')
    assert bf==FIXED and len(rows)==len(baseline)==1573
    assert [r['insref_id'] for r in rows]==[r['insref_id'] for r in baseline]
    assert len({r['insref_id'] for r in rows})==1573
    for r,b in zip(rows,baseline):
        for f in FIXED: assert r[f]==b[f],(r['insref_id'],f,r[f],b[f])
    spec=json.loads((BASE/'model_specification.json').read_text())
    columns=spec['columns'];n=len(rows);p=len(columns)
    X=np.array([[1.]+[float(int(r['annual_installments'])==a) for a in range(3,13) if a!=8]+[float(r['provisional_none'])]+[float(r[c]) for c in COV]+[float(r['prefecture_code']==f'{a:02d}') for a in range(2,48)] for r in rows])
    y=np.array([float(r['collection_ratio_fraction']) for r in rows])
    xsha=hashlib.sha256(X.astype('<f8').tobytes()).hexdigest();ysha=hashlib.sha256(y.astype('<f8').tobytes()).hexdigest()
    assert xsha==spec['X_memory_sha256'] and ysha==spec['y_memory_sha256']
    assert np.isfinite(X).all() and np.isfinite(y).all() and X.shape==(1573,61)
    N=np.array([float(r[args.count_column]) for r in rows]);A=np.array([float(r[args.ordinary_assessment_column]) for r in rows]);T=np.array([float(r[args.total_assessment_column]) for r in rows])
    assert all(np.isfinite(v).all() and (v>0).all() for v in [N,A,T])
    assert (N==np.round(N)).all(),'N not integer'
    assert (A<=T).all(),'Ordinary assessment exceeds total assessment: inspect definition'
    ln_error=float(np.max(abs(np.log(N)-X[:,columns.index('ln_insured_end_fy_count')])));assert ln_error<1e-12
    weights=[np.ones(n),A,N,N*A/T]
    annual=np.array([int(r['annual_installments']) for r in rows]);none=np.array([int(r['provisional_none']) for r in rows]);pref=np.array([r['prefecture_code'] for r in rows])
    groups=[('national','all',np.ones(n,dtype=bool))]
    groups += [('annual',str(a),annual==a) for a in range(3,13)]
    groups += [('classification','yes' if d==0 else 'no',none==d) for d in [0,1]]
    groups += [('cell',f'{a}_'+('yes' if d==0 else 'no'),(annual==a)&(none==d)) for a in range(3,13) for d in [0,1]]
    groups += [('prefecture',f'{a:02d}',pref==f'{a:02d}') for a in range(1,48)]
    group_results=[];national=[];coefs=[];contrasts=[];fits=[];rowdiag=[];rawcomp=[];quantiles=[];betas=[];definitions=[]
    _,fixedcoef=read(BASE/'fixed_ols_coefficients_all61.csv');fixed=np.array([float(r['estimate_fraction']) for r in fixedcoef]);assert [r['term'] for r in fixedcoef]==columns
    for name,w in zip(SCHEMES,weights):
        assert np.isfinite(w).all() and (w>0).all()
        wn=w/w.mean();shares=w/w.sum();root=np.sqrt(wn);Z=X*root[:,None];v=y*root
        beta,_,rank,s=np.linalg.lstsq(Z,v,rcond=None)
        bqr,_,rankqr,_=linalg.lstsq(Z,v,lapack_driver='gelsy',check_finite=True)
        assert rank==rankqr==p,(name,rank,rankqr)
        cross=float(np.max(abs(beta-bqr)));assert cross<1e-11,(name,cross)
        # Independent scale-invariance check at a different positive global normalization.
        bscale=np.linalg.lstsq(X*np.sqrt(shares)[:,None],y*np.sqrt(shares),rcond=None)[0]
        scale_err=float(np.max(abs(beta-bscale)));assert scale_err<1e-11
        if name=='equal_insurer':assert np.max(abs(beta-fixed))<1e-12
        betas.append(beta)
        fitted=X@beta;resid=y-fitted
        q,_=np.linalg.qr(Z,mode='reduced');h=np.sum(q*q,axis=1)
        assert abs(h.sum()-p)<1e-10 and h.min()>0 and h.max()<1+1e-10
        norm_resid=float(np.max(abs(X.T@(wn*resid))));assert norm_resid<1e-7
        desc=np.sort(shares)[::-1]
        nr={'scheme':name,'n':n,'positive_weights':int((w>0).sum()),'raw_weight_sum':scalar(w.sum()),'mean1_scaling_divisor':scalar(w.mean()),'weighted_mean_fraction':wmean(y,w),'weighted_mean_percent':100*wmean(y,w),'HHI':scalar(np.dot(shares,shares)),'Kish_effective_n':scalar(1/np.dot(shares,shares)),'Kish_n_over_n':scalar(1/np.dot(shares,shares)/n),'min_share':scalar(shares.min()),'max_share':scalar(shares.max())}
        for k in [1,5,10]: nr[f'top_{k}_share']=scalar(desc[:k].sum())
        for pct in [1,5,10]:
            k=int(np.ceil(n*pct/100));nr[f'top_{pct}pct_n']=k;nr[f'top_{pct}pct_share']=scalar(desc[:k].sum())
        national.append(nr)
        definitions.append({'scheme':name,'raw_formula':{'equal_insurer':'1','ordinary_current_assessment':'A','all_first_category_insured':'N','proxy_N_times_A_over_T':'N*A/T'}[name],'N_column':args.count_column,'A_column':args.ordinary_assessment_column,'T_column':args.total_assessment_column,'numerical_scaling':'raw w / arithmetic mean(raw w); all 1573 rows retained','target_interpretation':{'equal_insurer':'Equal insurer; baseline target','ordinary_current_assessment':'Ordinary current assessment monetary scale; not insurer or person target','all_first_category_insured':'All first-category insured year-end population scale; not ordinary-only people','proxy_N_times_A_over_T':'Constructed exposure proxy; NOT observed ordinary-insured count'}[name]})
        for u in [0,.01,.05,.25,.5,.75,.95,.99,1]:
            quantiles.append({'scheme':name,'quantile':u,'raw_weight':scalar(np.quantile(w,u)),'mean1_weight':scalar(np.quantile(wn,u)),'weight_share':scalar(np.quantile(shares,u)),'quantile_method':'numpy linear, unweighted insurer distribution'})
        for kind,label,mask in groups:
            k=int(mask.sum());wg=w[mask];total=float(wg.sum())
            item={'scheme':name,'group_type':kind,'group':label,'insurer_n':k,'raw_weight_total':total,'share_of_national_weight':total/w.sum(),'within_group_Kish_n':float(total**2/np.dot(wg,wg)) if k else None,'largest_within_group_weight_share':float(wg.max()/total) if k else None,'weighted_mean_fraction':wmean(y[mask],wg) if k else None,'weighted_mean_percent':100*wmean(y[mask],wg) if k else None,'empty_cell':k==0,'undefined_mean_reason':'empty group' if k==0 else ''}
            group_results.append(item)
        for j,term in enumerate(columns):
            coefs.append({'scheme':name,'column_1_based':j+1,'term':term,'estimate_fraction':scalar(beta[j]),'estimate_pp':scalar(100*beta[j]),'equal_insurer_estimate_pp':scalar(100*fixed[j]),'change_from_equal_pp':scalar(100*(beta[j]-fixed[j])),'sign_reversal_vs_equal':bool(beta[j]*fixed[j]<0),'inference':'point estimate only; no interval or test'})
        bnone=beta[columns.index('provisional_none')]
        ctr=[('provisional_no_minus_yes',bnone)]
        for a in range(3,13):
            ba=0. if a==8 else beta[columns.index(f'annual_{a}')]
            ctr.append((f'annual_{a}_minus_8_same_classification',ba))
            ctr.append((f'model_{a}yes_minus_8no',ba-bnone))
        for label,val in ctr:contrasts.append({'scheme':name,'contrast':label,'estimate_fraction':scalar(val),'estimate_pp':scalar(100*val),'inference':'additive model point contrast; no interval or test'})
        for label,ma,mb in [('no_minus_yes',none==1,none==0),('12yes_minus_8no',(annual==12)&(none==0),(annual==8)&(none==1))]:
            ga=wmean(y[ma],w[ma]);gb=wmean(y[mb],w[mb]);rawcomp.append({'scheme':name,'contrast':label,'group_A_n':int(ma.sum()),'group_B_n':int(mb.sum()),'group_A_mean_percent':100*ga,'group_B_mean_percent':100*gb,'raw_difference_pp':100*(ga-gb),'inference':'unadjusted weighted group means only'})
        fits.append({'scheme':name,'n':n,'p':p,'rank_SVD':int(rank),'rank_QR':int(rankqr),'condition_number_original_columns':scalar(s[0]/s[-1]),'minimum_singular_value_mean1_weight':scalar(s[-1]),'maximum_QR_SVD_difference_fraction':cross,'maximum_global_weight_rescaling_difference_fraction':scale_err,'maximum_abs_normal_equation_residual':norm_resid,'weighted_residual_sum':scalar(np.dot(wn,resid)),'weighted_R2':float(1-np.dot(wn,resid**2)/np.dot(wn,(y-wmean(y,wn))**2)),'weighted_RMSE_fraction':float(np.sqrt(np.dot(wn,resid**2)/wn.sum())),'fitted_min_fraction':scalar(fitted.min()),'fitted_max_fraction':scalar(fitted.max()),'fitted_below0_n':int((fitted<0).sum()),'fitted_above1_n':int((fitted>1).sum()),'leverage_sum':scalar(h.sum()),'leverage_max':scalar(h.max()),'leverage_over_2p_n':int((h>2*p/n).sum()),'leverage_over_3p_n':int((h>3*p/n).sum()),'leverage_ge_0_999999_n':int((h>=.999999).sum()),'new_inference':False})
        for i,r in enumerate(rows):rowdiag.append({'scheme':name,'insref_id':r['insref_id'],'prefecture_code':r['prefecture_code'],'annual_installments':r['annual_installments'],'provisional_none':r['provisional_none'],'raw_weight':scalar(w[i]),'mean1_weight':scalar(wn[i]),'national_weight_share':scalar(shares[i]),'observed_fraction':scalar(y[i]),'fitted_fraction':scalar(fitted[i]),'residual_fraction':scalar(resid[i]),'weighted_design_leverage':scalar(h[i])})
    B=np.stack(betas)
    ranges=[{'term':term,'minimum_estimate_pp':scalar(100*B[:,j].min()),'maximum_estimate_pp':scalar(100*B[:,j].max()),'range_pp':scalar(100*np.ptp(B[:,j])),'schemes_positive':int((B[:,j]>0).sum()),'schemes_negative':int((B[:,j]<0).sum()),'zero_crossing_across_weights':bool(B[:,j].min()<0<B[:,j].max())} for j,term in enumerate(columns)]
    scale_fields={'ordinary_assessment_per_all_insured_A_N':A/N,'ordinary_share_of_current_assessment_A_T':A/T,'total_current_assessment_per_all_insured_T_N':T/N,'all_first_category_insured_N':N}
    scale_desc=[]
    for kind,label,mask in groups:
        if kind=='prefecture':continue # Focused additional question only: policy groups/cells and nationwide.
        for field,z in scale_fields.items():
            k=int(mask.sum());zz=z[mask]
            for sname,ww in [('equal_insurer',np.ones(n)),('all_first_category_insured',N)]:
                scale_desc.append({'scheme':sname,'group_type':kind,'group':label,'field':field,'n':k,'weighted_mean':wmean(zz,ww[mask]) if k else None,'unweighted_min':scalar(zz.min()) if k else None,'unweighted_q25':scalar(np.quantile(zz,.25)) if k else None,'unweighted_median':scalar(np.median(zz)) if k else None,'unweighted_q75':scalar(np.quantile(zz,.75)) if k else None,'unweighted_max':scalar(zz.max()) if k else None,'percentiles_note':'Insurer-distribution percentiles repeated for transparency; mean uses named weight','units':'count for N; fraction for A/T; source monetary units per all-insured person for A/N,T/N'})
    U=np.log(A/N);V=np.log(A/T);W=np.log(T/N);L=np.log(N)
    decomp=[];corr=[]
    for name,w in [('equal_insurer',np.ones(n)),('all_first_category_insured',N)]:
        vu=wcov(U,U,w);vv=wcov(V,V,w);vw=wcov(W,W,w);cvw=wcov(V,W,w)
        decomp.append({'scheme':name,'variance_log_A_N':vu,'variance_log_A_T':vv,'variance_log_T_N':vw,'twice_covariance_components':2*cvw,'sum_components':vv+vw+2*cvw,'identity_error':vu-(vv+vw+2*cvw),'share_variance_log_A_T':vv/vu,'share_variance_log_T_N':vw/vu,'share_twice_covariance':2*cvw/vu,'interpretation':'algebraic descriptive decomposition, not causal or mediation'})
        for na,aa,nb,bb in [('log_A_N',U,'log_A_T',V),('log_A_N',U,'log_T_N',W),('log_A_T',V,'log_T_N',W),('log_N',L,'log_A_N',U)]:corr.append({'scheme':name,'variable_1':na,'variable_2':nb,'weighted_Pearson_correlation':wcorr(aa,bb,w),'inference':'descriptive only'})
    identity={'max_abs_log_identity_error':scalar(np.max(abs(U-(V+W)))),'assessment_minus_count_weighted_mean_pp':100*(wmean(y,A)-wmean(y,N)),'covariance_identity_rhs_pp':100*wcov(y,A/N,N)/wmean(A/N,N),'identity_error_pp':100*((wmean(y,A)-wmean(y,N))-wcov(y,A/N,N)/wmean(A/N,N)),'interpretation':'Exact reweighting identity; covariance is under N-weight distribution, no causal implication'}
    with localcontext() as ctx:
        ctx.prec=100
        dy=[Decimal(r['collection_ratio_fraction']) for r in rows];dn=[Decimal(r[args.count_column]) for r in rows];da=[Decimal(r[args.ordinary_assessment_column]) for r in rows];dt=[Decimal(r[args.total_assessment_column]) for r in rows]
        dws=[[Decimal(1)]*n,da,dn,[ni*ai/ti for ni,ai,ti in zip(dn,da,dt)]]
        decimal_results=[]
        for name,dw,nr in zip(SCHEMES,dws,national):
            den=sum(dw);num=sum((wi*yi for wi,yi in zip(dw,dy)),Decimal(0));mu=num/den
            decimal_results.append({'scheme':name,'exact_decimal_weight_sum':str(den),'weighted_rate_numerator_from_stored_y':str(num),'weighted_mean_fraction_decimal100':str(mu),'float_minus_decimal_fraction':str(Decimal(str(nr['weighted_mean_fraction']))-mu),'precision_note':'100-digit Decimal arithmetic over stored strings; proxy division and outcome source may already be rounded'})
        pooled={'ordinary_assessment_sum_exact':str(sum(da)),'all_first_category_insured_sum_exact':str(sum(dn)),'current_total_assessment_sum_exact':str(sum(dt)),'ordinary_collection_column':args.ordinary_collection_column,'status':'Numerator column unavailable; identity conditional on audited J/I outcome definition' if not args.ordinary_collection_column else 'Verified source numerator and denominator at stored precision'}
        if args.ordinary_collection_column:
            dj=[Decimal(r[args.ordinary_collection_column]) for r in rows]
            pooled['ordinary_collection_sum_exact']=str(sum(dj));pooled['pooled_J_over_A_fraction_decimal100']=str(sum(dj)/sum(da));pooled['assessment_weighted_mean_minus_pooled_fraction']=str(sum(ai*yi for ai,yi in zip(da,dy))/sum(da)-sum(dj)/sum(da));pooled['maximum_abs_row_y_minus_J_over_A']=str(max(abs(yi-ji/ai) for yi,ji,ai in zip(dy,dj,da)))
    tables=[('WEIGHT_DEFINITIONS4.csv',definitions),('NATIONAL_WEIGHT_SUMMARY4.csv',national),('WEIGHT_QUANTILES36.csv',quantiles),('GROUP_CELL_WEIGHT_DISTRIBUTIONS.csv',group_results),('ALL61_COEFFICIENTS_FOUR_WEIGHTS.csv',coefs),('FIXED_POLICY_CONTRASTS.csv',contrasts),('RAW_COMPARISONS8.csv',rawcomp),('FIT_AND_DESIGN_QA4.csv',fits),('ROW_WEIGHT_FIT_LEVERAGE.csv',rowdiag),('COEFFICIENT_RANGES61.csv',ranges),('SCALE_DESCRIPTIVES.csv',scale_desc),('LOG_SCALE_CORRELATIONS8.csv',corr),('LOG_SCALE_VARIANCE_DECOMPOSITION2.csv',decomp),('HIGH_PRECISION_WEIGHTED_MEANS4.csv',decimal_results)]
    for fn,rr in tables:save(fn,rr)
    dump('REWEIGHTING_IDENTITY.json',identity);dump('POOLED_COLLECTION_RATIO_CHECK.json',pooled)
    qa={'status':'PASS','started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,'python_executable':sys.executable,'platform':platform.platform(),'input_path':str(args.input),'input_sha256':sha(args.input),'baseline_csv_sha256':sha(BASE/'analysis_rows_1573.csv'),'model_specification_sha256':sha(BASE/'model_specification.json'),'baseline_string_fields_equal':True,'X_sha256':xsha,'y_sha256':ysha,'n':n,'p':p,'max_lnN_match_error':ln_error,'primary_solver':'numpy.linalg.lstsq SVD, rcond=None','solver_crosscheck':'scipy.linalg.lstsq pivoted QR, gelsy','weight_scaling':'divide raw weight by sample mean; rescaled-to-sum-1 independently cross-checked','rank_failures':[],'excluded_rows':0,'missing_or_nonpositive_weights':0,'imputation':False,'clipping':False,'new_inference':False,'all_tables':[{ 'file':fn,'rows':len(rr),'sha256':sha(HERE/fn)} for fn,rr in tables],'script_sha256':sha(__file__)}
    dump('RUN_RECEIPT.json',qa)
    print(json.dumps({'receipt':qa,'national':national,'raw_comparisons':rawcomp,'key_contrasts':[r for r in contrasts if r['contrast'] in ['provisional_no_minus_yes','annual_12_minus_8_same_classification','model_12yes_minus_8no']],'scale_decomposition':decomp,'scale_correlations':corr,'reweighting_identity':identity},indent=2,ensure_ascii=False))

if __name__=='__main__':
    try:main()
    except Exception as exc:
        dump('FAILED_ATTEMPT.json',{'utc':datetime.now(timezone.utc).isoformat(),'exception_type':type(exc).__name__,'message':str(exc),'arguments':sys.argv,'source_files_intentionally_modified':False})
        raise
