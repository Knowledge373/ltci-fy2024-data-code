#!/usr/bin/env python3
"""Compute fixed-input descriptives afresh using unchanged Decimal50 functions."""
from pathlib import Path
import argparse,csv,hashlib,json
from decimal import Decimal
from core import describe
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
 inp=ROOT/'data/common_input_1573.csv';assert hashlib.sha256(inp.read_bytes()).hexdigest()=='ea25bf07ed1ff42626ca73719ca252af7f563ca35cd19029139863cdbca468f2'
 with inp.open(encoding='utf-8-sig',newline='') as f:rr=list(csv.DictReader(f))
 ids=[r['insref_id'] for r in rr];assert len(ids)==len(set(ids))==1573 and ids==sorted(ids)
 v={r['insref_id']:{'rate':Decimal(r['collection_ratio_fraction']),'assessment':Decimal(r['ordinary_current_assessment_yen']),'collection':Decimal(r['ordinary_current_collected_yen'])} for r in rr}
 cells=[describe([r['insref_id'] for r in rr if int(r['annual_installments'])==a and r['provisional_none']==str(int(c=='なし'))],v,a,c) for a in range(3,13) for c in ['あり','なし']]
 for name,rows in [('DESCRIPTIVE_NATIONAL1',[describe(ids,v)]),('DESCRIPTIVE_CELLS20',cells)]:
  (out/(name+'.json')).write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
  with (out/(name+'.csv')).open('w',encoding='utf-8-sig',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 print('Computed national and 20 cell descriptives from fixed input.')
if __name__=='__main__':main()
