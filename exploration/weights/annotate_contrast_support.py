#!/usr/bin/env python3
"""Annotate fixed model contrasts with observed-cell support."""
import csv,json,hashlib,os
from pathlib import Path
from datetime import datetime,timezone
P=Path(os.environ['EXPLORATION_OUTPUT_DIR'])
I=Path(__file__).resolve().parents[1]/'inputs'/'common_input_1573.csv'
def read(f):return list(csv.DictReader(f.open(encoding='utf-8-sig')))
def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()
rows=read(I);out=[]
counts={(a,d):sum(int(r['annual_installments'])==a and int(r['provisional_none'])==d for r in rows) for a in range(3,13) for d in [0,1]}
for old in read(P/'FIXED_POLICY_CONTRASTS.csv'):
 r=dict(old);c=r['contrast'];target=None
 if c.startswith('model_'):
  a=int(c.split('_')[1].split('yes')[0]);target=counts[a,0]
  r.update(contrast_kind='additive_model_joint_cell_prediction',target_cell=f'{a}_yes',target_cell_observed_n=target,reference_cell='8_no',reference_cell_observed_n=counts[8,1],target_cell_absent=target==0,observed_same_classification_support='not applicable: classification differs',support_note=('Target cell is unobserved. This is an additive-model extrapolation, not an observed-cell comparison.' if target==0 else 'Both named cells are observed. This remains a fixed-covariate additive-model contrast, not the raw observed-cell mean difference.'))
 elif c.startswith('annual_'):
  a=int(c.split('_')[1]);support=[('yes' if d==0 else 'no') for d in [0,1] if counts[a,d] and counts[8,d]]
  r.update(contrast_kind='additive_model_annual_main_effect',target_cell=f'{a}_same_classification',target_cell_observed_n='',reference_cell='8_same_classification',reference_cell_observed_n='',target_cell_absent='',observed_same_classification_support=';'.join(support),support_note='Common annual coefficient under no interaction. Listed classifications have at least one observed row in both annual categories; this does not establish adequate overlap or a direct group-mean effect.')
 else:
  r.update(contrast_kind='additive_model_common_classification_effect',target_cell='no_at_common_annual_and_background',target_cell_observed_n='',reference_cell='yes_at_common_annual_and_background',reference_cell_observed_n='',target_cell_absent='',observed_same_classification_support='not applicable',support_note='Common no-minus-yes coefficient pooled across annual categories. Only 4,6,7,8,9,10 contain both classifications; annual categories 3,5,11,12 have one classification absent.')
 r['direct_observed_comparison_or_causal_effect']=False
 out.append(r)
f=P/'FIXED_POLICY_CONTRASTS_WITH_SUPPORT_ANNOTATION.csv'
with f.open('w',encoding='utf-8-sig',newline='') as fp:
 w=csv.DictWriter(fp,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
assert sum(r['target_cell_absent'] is True for r in out)==8
print(json.dumps({'status':'PASS','annotated_rows':len(out),'unobserved_target_cell_rows':8}))
