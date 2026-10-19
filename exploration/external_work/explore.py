#!/usr/bin/env python3
"""Fixed all-row and top-assessment-exclusion exploratory calculations."""
import argparse, csv, hashlib, io, json, math, platform, sys, zipfile
from collections import Counter
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from pathlib import Path
import numpy as np

WEIGHTS = ['equal', 'ordinary_assessment', 'all_insured_N', 'proxy_N_Ao_over_T']
POLICY = ['annual_4','annual_5','annual_6','annual_7','annual_9','annual_10','annual_12','provisional_none']
EXPECTED_CSV = 'ea25bf07ed1ff42626ca73719ca252af7f563ca35cd19029139863cdbca468f2'

def dump(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')

def save(path, rows):
    if not rows: return
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def digest(b): return hashlib.sha256(b).hexdigest()
def csvrows(b): return list(csv.DictReader(io.StringIO(b.decode('utf-8-sig'))))
def q(a, p): return float(np.quantile(a,p,method='linear'))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--input',type=Path,required=True);ap.add_argument('--model-specification',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args();out=args.out
    out.mkdir(parents=True,exist_ok=False)
    started=datetime.now(timezone.utc).isoformat()
    cb=args.input.read_bytes(); assert digest(cb)==EXPECTED_CSV
    rows=csvrows(cb);fields=list(rows[0]);n=len(rows);assert n==1573 and len(fields)==15
    ids=[r['insref_id'] for r in rows];assert len(set(ids))==n and ids==sorted(ids)
    missing={k:sum(not r[k].strip() or r[k].lower() in ['na','nan','null','none'] for r in rows) for k in fields}
    assert not any(missing.values())
    for r in rows:
        assert r['fiscal_year']=='2024' and r['provisional_none'] in ['0','1']
        assert 3<=int(r['annual_installments'])<=12 and 1<=int(r['prefecture_code'])<=47
        assert all(math.isfinite(float(r[k])) for k in fields[1:])
    spec=json.loads(args.model_specification.read_text())
    pref=np.array([r['prefecture_code'] for r in rows]); annual=np.array([int(r['annual_installments']) for r in rows]); d=np.array([int(r['provisional_none']) for r in rows])
    y=np.array([float(r['collection_ratio_fraction']) for r in rows]);assert ((0<=y)&(y<=1)).all()
    N=np.array([int(r['insured_count_end_of_fiscal_year']) for r in rows],dtype=np.int64)
    Ao=np.array([int(r['ordinary_current_assessment_yen']) for r in rows],dtype=np.int64)
    J=np.array([int(r['ordinary_current_collected_yen']) for r in rows],dtype=np.int64)
    As=np.array([int(r['special_current_assessment_yen']) for r in rows],dtype=np.int64)
    T=np.array([int(r['total_current_assessment_yen']) for r in rows],dtype=np.int64)
    assert (N>0).all() and (Ao>0).all() and (J>0).all() and (As>=0).all() and np.array_equal(T,Ao+As)
    X=np.column_stack([np.ones(n) if term=='intercept' else (annual==int(term.split('_')[1])).astype(float) if term.startswith('annual_') else
        (pref==term.split('_')[1]).astype(float) if term.startswith('prefecture_') else np.array([float(r[term]) for r in rows]) for term in spec['columns']])
    assert digest(X.astype('<f8').tobytes())==spec['X_memory_sha256'] and digest(y.astype('<f8').tobytes())==spec['y_memory_sha256']
    with localcontext() as c:
        c.prec=28
        assert all(Decimal(r['collection_ratio_fraction'])==Decimal(r['ordinary_current_collected_yen'])/Decimal(r['ordinary_current_assessment_yen']) for r in rows)
    with localcontext() as c:
        c.prec=80
        proxy=np.array([float(Decimal(int(nn))*Decimal(int(aa))/Decimal(int(tt))) for nn,aa,tt in zip(N,Ao,T)])
    weights={'equal':np.ones(n),'ordinary_assessment':Ao.astype(float),'all_insured_N':N.astype(float),'proxy_N_Ao_over_T':proxy}
    ordidx=sorted(range(n),key=lambda i:(-int(Ao[i]),ids[i]));top16=ordidx[:math.ceil(.01*n)]
    masks={'all':np.ones(n,dtype=bool),'exclude_Ao_top16':np.array([i not in set(top16) for i in range(n)])}
    bounds=np.quantile(N,[.25,.5,.75],method='linear');quart=np.searchsorted(bounds,N,side='left')+1
    audit={'n':n,'columns':fields,'missing_by_column':missing,'unique_ids':len(set(ids)),'prefectures':len(set(pref)),
           'special_zero_n':int((As==0).sum()),'observed_cells':len(set(zip(annual.tolist(),d.tolist()))),'empty_cells':20-len(set(zip(annual.tolist(),d.tolist()))),
           'all_y_match_J_over_Ao_Decimal28':True,'X_y_hash_match_spec':True,
           'original_source_reconstruction_performed':False,'N_quartile_bounds_persons':bounds.tolist(),
           'top16_rule':'Ao descending; ties ID ascending; same16 for all weights','top16_ids':[ids[i] for i in top16],
           'all15_fields_used':'IDs/year checks; prefecture and model inputs; N,Ao,As,T for weights/descriptions; J for accounting identity',
           'input_csv_sha256':EXPECTED_CSV}
    dump(out/'input_audit.json',audit)
    save(out/'fixed_exclusion16.csv',[{'insref_id':ids[i],'prefecture_code':pref[i],'annual_installments':int(annual[i]),'provisional_none':int(d[i]),'ordinary_current_assessment_yen':int(Ao[i])} for i in top16])
    trials=[];fits=[];policy=[];means=[];concentration=[];groupmeans=[];all_betas={}
    for target,mask in masks.items():
        for weight in WEIGHTS:
            w=weights[weight][mask];a=w/w.sum();xx=X[mask];yy=y[mask];sw=np.sqrt(w/w.mean())
            run=f'{target}__{weight}'
            try:
                beta,_,rank,s=np.linalg.lstsq(xx*sw[:,None],yy*sw,rcond=None)
                ok=rank==61 and np.isfinite(beta).all()
                status='computed' if ok else 'withheld_rank_or_nonfinite'
                trials.append({'trial_id':run,'kind':'WLS_point_only','target':target,'weight':weight,'n':int(mask.sum()),'columns':61,'rank':int(rank),'status':status,'new_inference':False,'failure':'' if ok else 'rank below61 or nonfinite coefficients'})
                if ok:
                    all_betas[run]=beta
                    for j,term in enumerate(spec['columns']):
                        rec={'trial_id':run,'target':target,'weight':weight,'term':term,'estimate_fraction':float(beta[j]),'estimate_pp':float(beta[j]*100)}
                        fits.append(rec)
                        if term in POLICY:policy.append(rec.copy())
                    pred=xx@beta
                    trials[-1].update({'prediction_below0':int((pred<0).sum()),'prediction_above1':int((pred>1).sum()),'unscaled_weighted_condition_number':float(s[0]/s[-1])})
            except np.linalg.LinAlgError as ex:
                trials.append({'trial_id':run,'kind':'WLS_point_only','target':target,'weight':weight,'n':int(mask.sum()),'columns':61,'rank':'','status':'failed','new_inference':False,'failure':str(ex)})
            means.append({'target':target,'weight':weight,'n':int(mask.sum()),'weight_sum':float(w.sum()),'weighted_mean_percent':float(a@yy*100)})
            aa=np.sort(a)[::-1];ps={str(p):float(a[pref[mask]==p].sum()) for p in sorted(set(pref[mask]))};pv=sorted(ps.values(),reverse=True)
            concentration.append({'target':target,'weight':weight,'n':int(mask.sum()),'top1_share_percent':float(aa[0]*100),'top16_share_percent':float(aa[:16].sum()*100),
                'inverse_sum_squared_normalized_weights':float(1/(a@a)),'prefecture_top1_share_percent':pv[0]*100,'prefecture_top5_share_percent':sum(pv[:5])*100,
                'special_zero_weight_share_percent':float(a[(As[mask]==0)].sum()*100),
                'Ao_excluded16_share_of_original_weight_percent':float(weights[weight][top16].sum()/weights[weight].sum()*100)})
            if target=='all':
                for dd in [0,1]:
                    gg=d==dd;wg=weights[weight][gg]
                    groupmeans.append({'weight':weight,'provisional_none':dd,'n':int(gg.sum()),'weighted_mean_percent':float(wg@y[gg]/wg.sum()*100)})
    save(out/'coefficients_all61.csv',fits);save(out/'policy8_points.csv',policy);save(out/'weighted_means.csv',means);save(out/'weight_concentration.csv',concentration);save(out/'classification_weighted_means.csv',groupmeans)
    deltas=[]
    for weight in WEIGHTS:
        k1='all__'+weight;k2='exclude_Ao_top16__'+weight
        if k1 in all_betas and k2 in all_betas:
            for term in POLICY:
                j=spec['columns'].index(term);deltas.append({'weight':weight,'term':term,'all_pp':float(all_betas[k1][j]*100),'exclude16_pp':float(all_betas[k2][j]*100),'exclude_minus_all_pp':float((all_betas[k2][j]-all_betas[k1][j])*100)})
    save(out/'exclusion_deltas_policy8.csv',deltas)
    dist=[]
    groups=[('national',np.ones(n,dtype=bool)),('provisional_yes',d==0),('provisional_no',d==1)]+[(f'N_quartile_{k}',quart==k) for k in range(1,5)]
    for label,m in groups:
        yy=y[m]*100
        rec={'group':label,'n':int(m.sum()),'mean_percent':float(yy.mean()),'median_percent':q(yy,.5),'q10_percent':q(yy,.1),'q25_percent':q(yy,.25),'q75_percent':q(yy,.75),'q90_percent':q(yy,.9),
             'below90_n':int((yy<90).sum()),'below90_share_percent':float((yy<90).mean()*100),'N_min':int(N[m].min()),'N_max':int(N[m].max()),'median_N':q(N[m],.5),
             'median_Ao_over_T_percent':q(Ao[m]/T[m],.5)*100,'median_T_per_N_yen':q(T[m]/N[m],.5)}
        dist.append(rec)
    save(out/'rate_distributions.csv',dist)
    cells=[]
    for count in range(3,13):
        for dd in [0,1]:
            m=(annual==count)&(d==dd);nn=int(m.sum())
            r={'annual_installments':count,'provisional_none':dd,'n':nn,'equal_mean_percent':float(y[m].mean()*100) if nn else '',
               'median_N':q(N[m],.5) if nn else '', 'median_Ao_over_T_percent':q(Ao[m]/T[m],.5)*100 if nn else '', 'median_T_per_N_yen':q(T[m]/N[m],.5) if nn else ''}
            for weight in WEIGHTS:
                w=weights[weight];r[weight+'_national_share_percent']=float(w[m].sum()/w.sum()*100)
                r[weight+'_within_cell_concentration_number']=float(w[m].sum()**2/(w[m]@w[m])) if nn else ''
            cells.append(r)
    save(out/'institutional_cells20.csv',cells)
    with localcontext() as c:
        c.prec=80
        sy=sum((Decimal(r['collection_ratio_fraction'])*Decimal(r['ordinary_current_assessment_yen']) for r in rows),Decimal(0))/sum(int(v) for v in Ao)
        rat=Decimal(sum(int(v) for v in J))/Decimal(sum(int(v) for v in Ao))
        identity={'Ao_weighted_saved_y_percent':str(sy*100),'sum_J_over_sum_Ao_percent':str(rat*100),'difference_pp':str((sy-rat)*100)}
    ap=weights['proxy_N_Ao_over_T']/weights['proxy_N_Ao_over_T'].sum();aa=Ao/Ao.sum()
    proxy_stats={'formula':'N*Ao/T = Ao/(T/N)','unit':'constructed population-equivalent proxy, not observed ordinary persons',
                 'normalized_proxy_vs_Ao_total_variation':float(.5*np.abs(ap-aa).sum()),
                 'Ao_over_T_quantiles_percent':{str(p):q(Ao/T,p)*100 for p in [.1,.25,.5,.75,.9]},
                 'T_over_N_quantiles_yen':{str(p):q(T/N,p) for p in [.1,.25,.5,.75,.9]},
                 'normalized_proxy_over_Ao_weight_ratio_quantiles':{str(p):q(ap/aa,p) for p in [.1,.25,.5,.75,.9]},'money_identity_check':identity}
    dump(out/'proxy_and_money_identity.json',proxy_stats)
    for e in ['E1_weight_comparison','E2_distribution_concentration','E3_institutional_cells']:
        trials.append({'trial_id':e,'kind':'descriptive','target':'all; declared exclusion sensitivity for weight summary','weight':'as declared','n':n,'columns':15,'rank':'','status':'computed','new_inference':False,'failure':''})
    dump(out/'trial_registry.json',trials)
    assert len([t for t in trials if t['kind']=='WLS_point_only'])==8
    assert digest(args.input.read_bytes())==EXPECTED_CSV
    receipt={'started_utc':started,'finished_utc':datetime.now(timezone.utc).isoformat(),'python':platform.python_version(),'numpy':np.__version__,'platform':platform.platform(),
             'method':'NumPy lstsq SVD, sqrt(weight/mean(weight)) scaling; type7 quantiles; no inference', 'fit_attempts':8,
             'successful_fits':sum(t['kind']=='WLS_point_only' and t['status']=='computed' for t in trials),'failed_or_withheld_fits':[t['trial_id'] for t in trials if t['status']!='computed'],
             'input_csv_unchanged':True,'new_uncertainty_estimates':False,'additional_interaction_models':False,'input_edits':False}
    dump(out/'run_receipt.json',receipt)
    print(json.dumps({'receipt':receipt,'input_audit':audit,'weighted_means':means,'policy8':policy,'distribution':dist,'proxy':proxy_stats,'concentration':concentration},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
