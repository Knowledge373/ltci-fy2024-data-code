#!/usr/bin/env python3
"""Fixed distribution and influence calculations for the supplied dataset."""
from __future__ import annotations
import os
for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
    os.environ[key]='1'
import argparse, csv, hashlib, json, platform, sys, traceback
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np
import pandas as pd
import scipy
from scipy import linalg, stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

POLICY=['annual_4','annual_5','annual_6','annual_7','annual_9','annual_10','annual_12','provisional_none']
COV=['income_low3_share','age_75_84_share','age_85plus_share','ln_insured_end_fy_count']
Q=[0,.01,.05,.10,.25,.50,.75,.90,.95,.99,1]
LEDGER=[]

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def utc(): return datetime.now(timezone.utc).isoformat()
def readcsv(p):
    with open(p,encoding='utf-8-sig',newline='') as f: return list(csv.DictReader(f))
def clean(obj):
    if isinstance(obj,dict): return {str(k):clean(v) for k,v in obj.items()}
    if isinstance(obj,(list,tuple,np.ndarray)): return [clean(v) for v in obj]
    if isinstance(obj,(np.integer,)): return int(obj)
    if isinstance(obj,(np.bool_,)): return bool(obj)
    if isinstance(obj,(float,np.floating)): return float(obj) if np.isfinite(obj) else None
    return obj
def dump(out,name,obj): (out/name).write_text(json.dumps(clean(obj),indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def save(out,name,rows): pd.DataFrame(rows).to_csv(out/name,index=False,encoding='utf-8-sig',float_format='%.17g')
def log(question,status,note): LEDGER.append({'time_utc':utc(),'question':question,'status':status,'note':note})
def groups(values,ids,k):
    order=np.lexsort((ids,values)); b=np.empty(len(values),int)
    b[order]=np.arange(len(values))*k//len(values)+1
    return b
def quant_summary(y):
    n=len(y)
    r={'n':n,'mean_percent':float(np.mean(y)*100) if n else None,
       'median_percent':float(np.median(y)*100) if n else None,
       'sd_pp':float(np.std(y,ddof=1)*100) if n>1 else None,
       'skewness_bias_corrected':float(stats.skew(y,bias=False)) if n>=3 and np.ptp(y)>0 else None,
       'dispersion_defined':n>1}
    for q in Q:r[f'quantile_{q:g}_percent']=float(np.quantile(y,q,method='linear')*100) if n else None
    for name,condition in [('lt90',y<.90),('lt95',y<.95),('ge99',y>=.99),('ge99p5',y>=.995),('eq100',y==1)]:
        r[name+'_n']=int(condition.sum());r[name+'_pct']=float(condition.mean()*100) if n else None
    return r

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--baseline',type=Path,required=True)
    ap.add_argument('--extended',type=Path,required=True)
    ap.add_argument('--extended-sha256',required=True)
    ap.add_argument('--assessed-column',required=True)
    ap.add_argument('--received-column',required=True)
    ap.add_argument('--population-column',required=True)
    ap.add_argument('--amount-unit',required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
    assert not (out/'COMPUTATION_STARTED.json').exists(),'Existing computation cannot be silently overwritten.'
    dump(out,'COMPUTATION_STARTED.json',{'started_utc':utc(),'command':sys.argv,'script_sha256':sha(__file__)})
    source_paths=[args.extended,args.baseline/'analysis_rows_1573.csv',args.baseline/'model_specification.json',
                  args.baseline/'fixed_ols_coefficients_all61.csv']
    manifest={str(p):sha(p) for p in source_paths}
    dump(out,'INPUT_HASHES_BEFORE.json',manifest)
    assert sha(args.extended)==args.extended_sha256,'Released extended hash mismatch'
    spec=json.loads((args.baseline/'model_specification.json').read_text())
    assert sha(args.baseline/'analysis_rows_1573.csv')==spec['analysis_csv_sha256']
    rows=readcsv(args.baseline/'analysis_rows_1573.csv');ext=readcsv(args.extended)
    ids=np.array([r['insref_id'] for r in rows]);n=len(ids)
    assert n==1573 and len(set(ids))==n and list(ids)==sorted(ids)
    assert len(ext)==n and len({r['insref_id'] for r in ext})==n
    er={r['insref_id']:r for r in ext};assert set(er)==set(ids)
    X=np.array([[1.]+[float(int(r['annual_installments'])==a) for a in range(3,13) if a!=8]+
                [float(r['provisional_none'])]+[float(r[c]) for c in COV]+
                [float(r['prefecture_code']==f'{p:02d}') for p in range(2,48)] for r in rows])
    y=np.array([float(r['collection_ratio_fraction']) for r in rows]);p=X.shape[1]
    xhash=hashlib.sha256(X.astype('<f8',copy=False).tobytes()).hexdigest()
    yhash=hashlib.sha256(y.astype('<f8',copy=False).tobytes()).hexdigest()
    assert xhash==spec['X_memory_sha256'] and yhash==spec['y_memory_sha256']
    assert X.shape==(n,61) and np.isfinite(X).all() and np.isfinite(y).all()
    assert ((y>=0)&(y<=1)).all()
    ann=np.array([int(r['annual_installments']) for r in rows]);prov=np.array([int(r['provisional_none']) for r in rows])
    ad=[Decimal(er[i][args.assessed_column]) for i in ids]
    rd=[Decimal(er[i][args.received_column]) for i in ids]
    popd=[Decimal(er[i][args.population_column]) for i in ids]
    assert all(v.is_finite() and v>0 for v in ad+popd)
    assert all(v.is_finite() and v>=0 for v in rd)
    assert all(r<=a for r,a in zip(rd,ad))
    assert all(v==v.to_integral_value() for v in popd),'Population must be a count'
    assessed=np.array([float(v) for v in ad]);received=np.array([float(v) for v in rd]);population=np.array([float(v) for v in popd])
    money_recon=np.array([float(r/a) for r,a in zip(rd,ad)])
    assert np.max(abs(money_recon-y))<1e-12,'Audited amounts do not reconstruct fixed outcome'
    assert np.max(abs(np.log(population)-X[:,14]))<1e-12,'Population differs from fixed ln count'
    b,_,rank,_=linalg.lstsq(X,y,lapack_driver='gelsy');bs,_,ranks,s=np.linalg.lstsq(X,y,rcond=None)
    assert rank==ranks==61
    ref=readcsv(args.baseline/'fixed_ols_coefficients_all61.csv')
    assert [r['term'] for r in ref]==spec['columns']
    bref=np.array([float(r['estimate_fraction']) for r in ref])
    assert np.max(abs(b-bref))<1e-12 and np.max(abs(b-bs))<1e-12
    fitted=X@b;e=y-fitted;mse=float(e@e/(n-p))
    qq,rr=np.linalg.qr(X,mode='reduced');w=linalg.solve_triangular(rr,qq.T,lower=False).T
    h=(qq**2).sum(axis=1);regular=h<1-1e-8
    assert h.min()>=0 and h.max()<1+1e-8
    student=np.full(n,np.nan);cook=np.full(n,np.nan)
    student[regular]=e[regular]/np.sqrt(mse*(1-h[regular]))
    cook[regular]=e[regular]**2*h[regular]/(p*mse*(1-h[regular])**2)
    js=[spec['columns'].index(t) for t in POLICY]
    delta=np.full((n,len(POLICY)),np.nan)
    delta[regular]=-100*e[regular,None]*w[regular][:,js]/(1-h[regular,None])
    pl=w[:,js]**2/(w[:,js]**2).sum(axis=0)
    assert np.max(abs(pl.sum(axis=0)-1))<1e-12
    fwlqa=[]
    for j,term in zip(js,POLICY):
        z=np.delete(X,j,axis=1);g,_,rnk,_=linalg.lstsq(z,X[:,j],lapack_driver='gelsy');assert rnk==60
        r=X[:,j]-z@g;fw=r/(r@r)
        fwlqa.append({'term':term,'weight_max_absolute_difference':float(np.max(abs(fw-w[:,j]))),
                      'estimate_absolute_difference_pp':float(abs(fw@y-b[j])*100)})
    assert max(r['weight_max_absolute_difference'] for r in fwlqa)<1e-10
    assert max(r['estimate_absolute_difference_pp'] for r in fwlqa)<1e-8
    save(out,'FWL_QA.csv',fwlqa)
    dump(out,'MODEL_QA.json',{'n':n,'p':p,'rank_QR':rank,'rank_SVD':ranks,'matrix_hash':xhash,'outcome_hash':yhash,
         'max_abs_fixed_coefficient_difference_fraction':np.max(abs(b-bref)),
         'max_abs_QR_SVD_coefficient_difference_fraction':np.max(abs(b-bs)),
         'max_abs_outcome_amount_reconstruction_difference_fraction':np.max(abs(money_recon-y)),
         'max_abs_log_count_difference':np.max(abs(np.log(population)-X[:,14])),
         'fitted_below_zero_n':(fitted<0).sum(),'fitted_above_one_n':(fitted>1).sum(),
         'unit_leverage_n':(~regular).sum(),'unit_leverage_ids':ids[~regular],
         'h_sum':h.sum(),'SSE_fraction_squared':e@e,'MSE_fraction_squared':mse,
         'source_rows_modified':0,'substantive_altered_sample_models':0,'new_inference':False})
    log('integrity','PASS','Audited shared hash, row identity, exact fixed X/y hashes, amounts/count reconstruction and all 61 coefficients pass.')

    # Q1: all predeclared group summaries and contrasts, including small/null results.
    summaries=[];sm={}
    masks=[('all',np.ones(n,bool))]+[(f'provisional_{d}',prov==d) for d in (0,1)]
    masks += [(f'annual_{a}',ann==a) for a in range(3,13)]
    masks += [(f'annual_{a}_provisional_{d}',(ann==a)&(prov==d)) for a in range(3,13) for d in (0,1)]
    for group,mask in masks:
        r={'group':group,**quant_summary(y[mask])};summaries.append(r);sm[group]=r
    save(out,'DISTRIBUTION_SUMMARY_ALL_GROUPS.csv',summaries)
    comparisons=[('no_minus_yes','provisional_1','provisional_0')]+[(f'{a}_minus_8',f'annual_{a}','annual_8') for a in range(3,13) if a!=8]
    comparisons += [('12_yes_minus_8_no','annual_12_provisional_0','annual_8_provisional_1')]
    contrasts=[]
    for label,a,c in comparisons:
        ma,mc=sm[a],sm[c]
        diff={key:ma[key]-mc[key] for key in ['mean_percent','median_percent','quantile_0.1_percent']}
        md=diff['mean_percent'];med=diff['median_percent'];q10=diff['quantile_0.1_percent']
        contrasts.append({'comparison':label,'group_A':a,'group_B':c,'n_A':ma['n'],'n_B':mc['n'],
                          'mean_difference_pp':md,'median_difference_pp':med,'q10_difference_pp':q10,
                          'median_minus_mean_contrast_pp':med-md,'q10_minus_mean_contrast_pp':q10-md,
                          'opposite_mean_median_sign':md*med<0,'prespecified_descriptive_screen':md*med<0 or abs(med-md)>=.5 or abs(q10-md)>=1})
    save(out,'DISTRIBUTION_CONTRASTS.csv',contrasts)
    log('Q1','COMPLETE',f'{len(summaries)} groups and {len(contrasts)} contrasts retained; empty and singleton cells explicit; no tests.')

    # Exact monetary sums retained; scale arrays used solely for bins/ranks and descriptive statistics.
    with localcontext() as ctx:
        ctx.prec=160;total_ad=sum(ad);total_remainder=sum(a-r for a,r in zip(ad,rd))
    def monetary(mask):
        ix=np.flatnonzero(mask)
        with localcontext() as ctx:
            ctx.prec=160;da=sum((ad[i] for i in ix),Decimal(0));ra=sum((rd[i] for i in ix),Decimal(0));remaining=da-ra
            return {'assessed_sum_exact':str(da),'received_sum_exact':str(ra),'remainder_sum_exact':str(remaining),
                    'pooled_monetary_collection_percent':float(ra/da*100) if da else None,
                    'share_total_assessed_pct':float(da/total_ad*100),
                    'share_total_remainder_pct':float(remaining/total_remainder*100) if total_remainder else None,
                    'amount_unit':args.amount_unit}
    scales={'assessed_amount':assessed,'insured_population':population}
    bins={};scale_rows=[];scale_screen=[];correlations=[]
    for name,values in scales.items():
        for k in (4,10):
            g=groups(values,ids,k);bins[(name,k)]=g
            for q in range(1,k+1):
                mask=g==q
                scale_rows.append({'scale':name,'bin_count':k,'bin':q,'scale_min':values[mask].min(),'scale_max':values[mask].max(),
                                   'scale_median':np.median(values[mask]),'mean_abs_residual_pp':np.mean(abs(e[mask]))*100,
                                   **quant_summary(y[mask]),**monetary(mask)})
        low=[r for r in scale_rows if r['scale']==name and r['bin_count']==10 and r['bin']==1][0]
        high=[r for r in scale_rows if r['scale']==name and r['bin_count']==10 and r['bin']==10][0]
        md=high['mean_percent']-low['mean_percent'];lt90=high['lt90_pct']-low['lt90_pct'];lt95=high['lt95_pct']-low['lt95_pct']
        scale_screen.append({'scale':name,'largest_minus_smallest_decile_mean_pp':md,
                             'largest_minus_smallest_decile_lt90_prevalence_pp':lt90,
                             'largest_minus_smallest_decile_lt95_prevalence_pp':lt95,
                             'prespecified_descriptive_screen':abs(md)>=1 or abs(lt90)>=5 or abs(lt95)>=5})
        for outcome,v in [('rate',y),('absolute_fixed_model_residual',abs(e))]:
            correlations.append({'scale':name,'outcome':outcome,'spearman_rho':float(stats.spearmanr(np.log(values),v).statistic),
                                 'p_value_computed_or_reported':False,'interpretation':'descriptive rank association; not causal'})
    save(out,'SCALE_BIN_SUMMARIES.csv',scale_rows);save(out,'SCALE_ENDPOINT_SCREENS.csv',scale_screen);save(out,'SCALE_RANK_CORRELATIONS.csv',correlations)
    cross=[]
    for a in range(1,5):
        for c in range(1,5):
            mask=(bins[('assessed_amount',4)]==a)&(bins[('insured_population',4)]==c)
            cross.append({'assessed_quartile':a,'population_quartile':c,**quant_summary(y[mask]),**monetary(mask)})
    save(out,'FISCAL_POPULATION_QUARTILE_CROSS.csv',cross)
    lowdec=groups(y,ids,10)==1
    save(out,'MONETARY_CONCENTRATION.csv',[{'group':'all',**monetary(np.ones(n,bool))},
        {'group':'lowest_outcome_decile',**monetary(lowdec)}]+[
        {'group':f'{name}_quartile_{q}',**monetary(bins[(name,4)]==q)} for name in scales for q in range(1,5)])
    log('Q2','COMPLETE','Both scales, all 28 scale bins, all 16 crossed cells, exact monetary sums and all rank correlations retained.')

    # Q3: predeclared concepts remain separate. Structural h=1 observations stay in data.
    maxdelta=np.full(n,np.nan);maxdelta[regular]=np.max(abs(delta[regular]),axis=1)
    metrics={'low_y':y,'absolute_residual':-abs(e),'whole_leverage':-h,'max_policy8_influence':-maxdelta}
    k=int(np.ceil(n*.01));sets={};rank_cols={}
    for name,v in metrics.items():
        order=np.lexsort((ids,np.where(np.isfinite(v),v,np.inf)));eligible=order[np.isfinite(v[order])]
        sets[name]=set(eligible[:k]);ranks=np.full(n,np.nan);ranks[eligible]=np.arange(len(eligible))+1;rank_cols[name]=ranks
    rowdiag=[]
    for i in range(n):
        r={key:rows[i][key] for key in ['insref_id','prefecture_code','annual_installments','provisional_none']}
        r.update({'outcome_percent':y[i]*100,'fitted_percent':fitted[i]*100,'residual_pp':e[i]*100,'whole_design_h':h[i],
                  'regular_deletion_formula_valid':regular[i],'undefined_reason':None if regular[i] else 'h>=1-1e-8; fixed-design row deletion rank deficient',
                  'internally_standardized_residual':student[i],'cook_distance':cook[i],
                  'max_absolute_policy8_deletion_delta_pp':maxdelta[i],
                  'assessed_amount':assessed[i],'received_amount':received[i],'amount_unit':args.amount_unit,'insured_population':population[i],
                  'flag_abs_internal_residual_gt3':bool(abs(student[i])>3) if regular[i] else None,
                  'flag_h_gt_2p_over_n':bool(h[i]>2*p/n),'flag_cook_gt_4_over_n_minus_p':bool(cook[i]>4/(n-p)) if regular[i] else None})
        for name in metrics:r['rank_'+name]=rank_cols[name][i];r['top1pct_'+name]=i in sets[name]
        for name in scales:
            for nb in (4,10):r[f'{name}_{nb}_bin']=bins[(name,nb)][i]
        for j,t in enumerate(POLICY):r[t+'_deletion_delta_pp']=delta[i,j];r[t+'_partial_leverage']=pl[i,j]
        rowdiag.append(r)
    save(out,'ROW_DIAGNOSTICS_1573.csv',rowdiag)
    save(out,'TOP1PCT_SET_MEMBERS.csv',[{'set':name,**rowdiag[i]} for name,sset in sets.items() for i in sorted(sset)])
    overlaps=[]
    for a,sa in sets.items():
        for c,sc in sets.items():
            if list(sets).index(a)<list(sets).index(c):
                overlaps.append({'set_A':a,'set_B':c,'set_A_n':len(sa),'set_B_n':len(sc),'intersection_n':len(sa&sc),
                                 'jaccard':len(sa&sc)/len(sa|sc),'intersection_ids':';'.join(ids[i] for i in sorted(sa&sc))})
    save(out,'TOP1PCT_SET_OVERLAPS.csv',overlaps)
    setsummary=[]
    for name,sset in sets.items():
        ix=np.array(sorted(sset));setsummary.append({'set':name,'n':len(ix),'median_assessed_amount':np.median(assessed[ix]),
            'median_insured_population':np.median(population[ix]),'median_outcome_percent':np.median(y[ix])*100,
            'median_absolute_residual_pp':np.median(abs(e[ix]))*100,'unit_leverage_n':int((~regular[ix]).sum())})
    save(out,'TOP1PCT_SET_SCALE_SUMMARY.csv',setsummary)
    flags=[r for r in rowdiag if r['flag_abs_internal_residual_gt3'] or r['flag_h_gt_2p_over_n'] or r['flag_cook_gt_4_over_n_minus_p']]
    save(out,'ALL_SCREENED_ROWS.csv',flags)
    concentration=[]
    masks2=[('lowest_outcome_decile',lowdec)]+[(f'{name}_quartile_{q}',bins[(name,4)]==q) for name in scales for q in range(1,5)]
    for j,term in enumerate(POLICY):
        absall=float(np.nansum(abs(delta[:,j])))
        for name,mask in masks2:
            concentration.append({'term':term,'group':name,'n':int(mask.sum()),'valid_deletion_n':int((mask&regular).sum()),
                                  'absolute_deletion_influence_sum_pp':float(np.nansum(abs(delta[mask,j]))),
                                  'share_total_absolute_deletion_influence_pct':float(np.nansum(abs(delta[mask,j]))/absall*100),
                                  'share_partial_leverage_pct':float(pl[mask,j].sum()*100)})
    save(out,'COEFFICIENT_INFLUENCE_CONCENTRATION.csv',concentration)
    coefsummary=[]
    for j,t in enumerate(POLICY):
        i=int(np.nanargmax(abs(delta[:,j])))
        coefsummary.append({'term':t,'fixed_estimate_pp':b[js[j]]*100,'maximum_abs_single_row_delta_pp':abs(delta[i,j]),
                            'maximum_delta_row':ids[i],'signed_maximum_delta_pp':delta[i,j],
                            'largest16_absolute_delta_share_pct':np.nansum(np.sort(abs(delta[regular,j]))[-k:])/np.nansum(abs(delta[:,j]))*100,
                            'undefined_unit_leverage_rows':int((~regular).sum())})
    save(out,'POLICY8_DIAGNOSTIC_SUMMARY.csv',coefsummary)
    selected={}
    for label,values in [('abs_residual',abs(e)),('h',h),('cook',cook),('annual_12_delta',abs(delta[:,POLICY.index('annual_12')])),('provisional_none_delta',abs(delta[:,POLICY.index('provisional_none')]))]:
        eligible=np.flatnonzero(regular);i=int(eligible[np.nanargmax(values[eligible])]);selected.setdefault(i,[]).append(label)
    qa=[]
    for i,reasons in sorted(selected.items()):
        mask=np.arange(n)!=i;bb,_,rank_del,_=linalg.lstsq(X[mask],y[mask],lapack_driver='gelsy')
        difference=(bb[js]-b[js])*100;err=float(np.max(abs(difference-delta[i])))
        qa.append({'insref_id':ids[i],'selection_reasons':';'.join(reasons),'n':int(mask.sum()),'columns_retained':p,
                   'rank':rank_del,'maximum_absolute_formula_refit_difference_pp':err,'passed':rank_del==p and err<1e-8})
        assert rank_del==p and err<1e-8
    save(out,'DIRECT_REFIT_FORMULA_QA.csv',qa)
    dump(out,'DIAGNOSTIC_COUNTS.json',{'top_set_n':k,'unit_leverage_n':int((~regular).sum()),
         'internal_standardized_residual_abs_gt3_n':int(np.nansum(abs(student)>3)),
         'h_gt_2p_over_n_n':int((h>2*p/n).sum()),'cook_gt_4_over_n_minus_p_n':int(np.nansum(cook>4/(n-p))),
         'h_screen_threshold':2*p/n,'cook_screen_threshold':4/(n-p),'formula_QA_refit_n':len(qa),
         'singular_deletion_handling':'Rows retained, ordinary standardized residual/Cook/deletion deltas undefined; no dropped singleton columns.'})
    log('Q3','COMPLETE',f'All {n} rows retained, {int((~regular).sum())} structural h=1 rows explicitly undefined for standard deletion formulas; {len(qa)} full-column direct-refit QA checks passed.')

    # Static scientific figure, purely descriptive and not a manuscript insertion.
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    histogram_lower=min(65,int(np.floor(y.min()*100)))
    axes[0,0].hist(y*100,bins=np.arange(histogram_lower,101),color='#397291',edgecolor='white')
    axes[0,0].set(xlabel='Collection rate (%)',ylabel='Insurers',title='A. Observed distribution (all 1,573 insurers)')
    for name,col,label in [('assessed_amount','#397291','Assessed amount'),('insured_population','#bc5a45','Insured population')]:
        selected_bins=[r for r in scale_rows if r['scale']==name and r['bin_count']==10]
        axes[0,1].plot(range(1,11),[r['mean_percent'] for r in selected_bins],marker='o',color=col,label=label)
        axes[1,0].plot(range(1,11),[r['lt95_pct'] for r in selected_bins],marker='o',color=col,label=label)
    axes[0,1].set(xlabel='Equal-count scale decile (small to large)',ylabel='Equal-insurer mean rate (%)',title='B. Mean rate by scale');axes[0,1].legend(frameon=False)
    axes[1,0].set(xlabel='Equal-count scale decile (small to large)',ylabel='Insurers with rate <95% (%)',title='C. Lower-rate prevalence by scale')
    show=regular;axes[1,1].scatter(h[show],abs(e[show])*100,c=y[show]*100,s=11,alpha=.55,cmap='viridis')
    axes[1,1].set(xlabel='Whole-design leverage h',ylabel='Absolute fixed-model residual (pp)',title='D. Leverage and residual surprise differ')
    axes[1,1].text(.02,.97,f'{int((~regular).sum())} structural h=1 singleton rows omitted\nfrom this panel only; retained in all analyses.',transform=axes[1,1].transAxes,va='top',fontsize=8)
    for ax in axes.flat:ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.15)
    fig.suptitle('Exploratory description, FY2024 • fixed 1,573-insurer dataset',fontsize=14)
    fig.savefig(out/'EXPLORATORY_DISTRIBUTION_SCALE_INFLUENCE.png',dpi=180)
    fig.savefig(out/'EXPLORATORY_DISTRIBUTION_SCALE_INFLUENCE.pdf');plt.close(fig)
    log('figure','COMPLETE','Four static descriptive panels; no confidence intervals; two structural leverage rows excluded only from scatter panel, labeled.')
    after={str(p):sha(p) for p in source_paths};assert after==manifest
    dump(out,'INPUT_HASHES_AFTER.json',after)
    dump(out,'SOFTWARE.json',{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,
        'scipy':scipy.__version__,'matplotlib':matplotlib.__version__,'platform':platform.platform(),'runtime_executable':sys.executable,
        'blas_threads':{key:os.environ[key] for key in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']}})
    save(out,'TRIAL_LEDGER.csv',LEDGER)
    dump(out,'COMPUTATION_COMPLETE.json',{'completed_utc':utc(),'status':'PASS','new_inference':False,'altered_sample_substantive_models':0,
         'input_hashes_unchanged':True,'script_sha256':sha(__file__),'output_files_n_before_manifest':len(list(out.iterdir()))})
    dump(out,'OUTPUT_HASHES.json',{p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name!='OUTPUT_HASHES.json'})
    print(json.dumps({'status':'PASS','output_dir':str(out),'overall':sm['all'],'screening':scale_screen,'overlaps':overlaps},ensure_ascii=False,indent=2))

if __name__=='__main__':
    try:main()
    except Exception as exc:
        print(traceback.format_exc(),file=sys.stderr)
        try:
            out=Path(sys.argv[sys.argv.index('--out')+1]);out.mkdir(parents=True,exist_ok=True)
            log('execution','FAILED',repr(exc));save(out,'TRIAL_LEDGER.csv',LEDGER)
            dump(out,'FAILURE.json',{'utc':utc(),'exception':repr(exc),'traceback':traceback.format_exc()})
        except Exception:pass
        raise
