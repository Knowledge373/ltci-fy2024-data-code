"""Original Decimal50 descriptive functions, extracted without numerical changes."""
from decimal import Decimal, localcontext
RATE_STATS = ['mean','min','q25','median','q75','max']

def quantile7(a, q):
    h=Decimal(len(a)-1)*q
    lo=int(h); weight=h-Decimal(lo)
    return a[lo] if lo==len(a)-1 else a[lo]+weight*(a[lo+1]-a[lo])

def describe(ids, values, annual=None, classification=None):
    out={'annual_installments':annual,'adopted_classification':classification,
         'n':len(ids),'n_valid_rate':len(ids),'n_valid_amounts':len(ids)}
    with localcontext() as c:
        c.prec=50
        if ids:
            ys=sorted(values[i]['rate'] for i in ids)
            stats=[sum(ys,Decimal(0))/Decimal(len(ys)),ys[0],quantile7(ys,Decimal('0.25')),
                   quantile7(ys,Decimal('0.5')),quantile7(ys,Decimal('0.75')),ys[-1]]
            for key,val in zip(RATE_STATS,stats):
                out[key+'_fraction']=str(val); out[key+'_percent']=str(val*100)
            assessment=sum((values[i]['assessment'] for i in ids),Decimal(0))
            collection=sum((values[i]['collection'] for i in ids),Decimal(0))
            out.update(sum_collection_yen=str(collection),sum_assessment_yen=str(assessment),
                       amount_aggregate_rate_fraction=str(collection/assessment),
                       amount_aggregate_rate_percent=str((collection/assessment)*100))
        else:
            for key in RATE_STATS: out[key+'_fraction']=out[key+'_percent']=None
            out.update(sum_collection_yen=None,sum_assessment_yen=None,
                       amount_aggregate_rate_fraction=None,amount_aggregate_rate_percent=None)
    return out
