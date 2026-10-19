#!/usr/bin/env python3
"""Recompute the fixed equal-insurer national model and two existing contrasts."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import argparse,csv,hashlib,json,platform
from decimal import Decimal,ROUND_HALF_EVEN
import numpy as np
import scipy
from scipy import linalg
from core import calc_model,contrast
ROOT=Path(__file__).resolve().parents[1]
TARGET=['annual_4','annual_5','annual_6','annual_7','annual_9','annual_10','annual_12','provisional_none','income_low3_share','age_75_84_share','age_85plus_share','ln_insured_end_fy_count']
def save(out,name,rr):
 (out/(name+'.json')).write_text(json.dumps(rr,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
 with (out/(name+'.csv')).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
 assert (platform.python_version(),np.__version__,scipy.__version__)==('3.12.14','2.3.5','1.17.0')
 inp=ROOT/'data/common_input_1573.csv';assert hashlib.sha256(inp.read_bytes()).hexdigest()=='ea25bf07ed1ff42626ca73719ca252af7f563ca35cd19029139863cdbca468f2'
 with inp.open(encoding='utf-8-sig',newline='') as f:rr=list(csv.DictReader(f))
 spec=json.loads((ROOT/'data/baseline_model_specification.json').read_text());cols=spec['columns'];cov=spec['covariates']
 X=np.asarray([[1.]+[float(int(r['annual_installments'])==n) for n in range(3,13) if n!=8]+[float(r['provisional_none'])]+[float(r[k]) for k in cov]+[float(int(r['prefecture_code'])==n) for n in range(2,48)] for r in rr],dtype='<f8')
 y=np.array([float(r['collection_ratio_fraction']) for r in rr],dtype='<f8');pref=np.array([int(r['prefecture_code']) for r in rr],dtype='<i4')
 assert X.shape==(1573,61) and hashlib.sha256(X.tobytes()).hexdigest()==spec['X_memory_sha256'] and hashlib.sha256(y.tobytes()).hexdigest()==spec['y_memory_sha256']
 idx=[cols.index(t) for t in TARGET];T=np.column_stack([(pref==g).astype(float) for g in range(1,48)]+[X[:,cols.index('annual_3')],X[:,cols.index('annual_11')]])
 Qt,_=linalg.qr(T,mode='economic');Rs=X[:,idx]-Qt@(Qt.T@X[:,idx]);ys=y-Qt@(Qt.T@y)
 full=calc_model(X,y,pref,cols,'full61');part=calc_model(Rs,ys,pref,TARGET,'partial12')
 assert np.allclose(full['beta'][idx],part['beta'],atol=1e-10,rtol=1e-8)
 assert np.allclose(full['V'][np.ix_(idx,idx)],part['V'],atol=1e-10,rtol=1e-8)
 assert np.allclose(full['df'][idx],part['df'],atol=1e-8,rtol=1e-8)
 save(out,'TARGET12_RESULTS',[{'term':t,**{f:float(full[k][i]) for f,k in [('estimate_fraction','beta'),('SE_fraction','se'),('df','df'),('CI_lower_fraction','lo'),('CI_upper_fraction','hi')]}} for t,i in zip(TARGET,idx)])
 save(out,'ALL61_COEFFICIENTS',[{'term':t,'estimate_fraction':float(b)} for t,b in zip(cols,full['beta'])])
 result=[];savedpoints={'10no_minus_12yes':Decimal('-0.0555453337598402'),'9no_minus_12yes':Decimal('-0.1626420458358717')}
 for annual in [10,9]:
  name=f'{annual}no_minus_12yes';c=np.zeros(61);c[cols.index(f'annual_{annual}')]=1;c[cols.index('provisional_none')]=1;c[cols.index('annual_12')]=-1
  z=contrast(full,c);q=contrast(part,c[idx]);assert np.allclose(z['P'],z['Pdirect'],atol=1e-10,rtol=1e-8)
  for field in ['estimate_fraction','variance_fraction_squared','SE_fraction','df','CI_lower_fraction','CI_upper_fraction']:assert np.isclose(z[field],q[field],atol=1e-8 if field=='df' else 1e-10,rtol=1e-8)
  pointpp=savedpoints[name];center=float(pointpp/Decimal(100));assert np.isclose(center,z['estimate_fraction'],atol=1e-10,rtol=1e-8)
  half=z['t975']*z['SE_fraction'];lo=center-half;hi=center+half
  result.append({'contrast':name,'saved_estimate_pp':str(pointpp),'estimate_fraction':center,'recomputed_estimate_fraction':z['estimate_fraction'],'variance_fraction_squared':z['variance_fraction_squared'],'SE_fraction':z['SE_fraction'],'SE_pp':z['SE_fraction']*100,'df':z['df'],'t975':z['t975'],'CI_lower_fraction':lo,'CI_upper_fraction':hi,'CI_lower_pp':lo*100,'CI_upper_pp':hi*100})
 save(out,'TWO_CONTRAST_RESULTS',result)
 money=[]
 for r in result:
  for Ao in [16000000,50000000,130000000]:
   point=Decimal(Ao)*Decimal(r['saved_estimate_pp'])/Decimal(100);low=Decimal(Ao)*Decimal(repr(r['CI_lower_pp']))/Decimal(100);high=Decimal(Ao)*Decimal(repr(r['CI_upper_pp']))/Decimal(100)
   def man(v):return str((v/Decimal(10000)).quantize(Decimal('.1'),rounding=ROUND_HALF_EVEN))
   money.append({'contrast':r['contrast'],'fixed_assessment_yen':Ao,'point_yen':str(point),'lower_yen':str(low),'upper_yen':str(high),'point_10k_yen':man(point),'lower_10k_yen':man(low),'upper_10k_yen':man(high)})
 save(out,'ILLUSTRATIVE_MONEY_SCALES',money)
 print('Fresh national CR2 and two contrasts; historical decimal reporting centers retained after fresh-fit agreement.')
if __name__=='__main__':main()
