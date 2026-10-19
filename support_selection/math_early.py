"""Fixed support-selection mathematical functions."""
from decimal import Decimal
from fractions import Fraction

VARIABLES = ['income_low3_share', 'age_75_84_share', 'age_85plus_share', 'ln_insured_end_fy_count']

CALIPERS = [('1', Fraction(1)), ('0.5', Fraction(1, 2))]

REASONS = ['missing_invalid_input', 'no_opposite_group_in_prefecture', 'no_joint_neighbor_within_fixed_caliper']

PRECISION = 80

def exact(x):
    return f'{x.numerator}/{x.denominator}'

def decimal_fraction(x):
    return Decimal(x.numerator) / Decimal(x.denominator)

def allowed_count(text):
    try:
        v = Fraction(text) if text not in (None, '') else None
    except (ValueError, TypeError, ZeroDivisionError):
        return None
    return v if v is not None and v.denominator == 1 and v >= 0 else None

def pooled_variance(a, b):
    if len(a) < 2 or len(b) < 2:
        return None
    z = a[0] * 0
    mean_a, mean_b = sum(a,z)/len(a), sum(b,z)/len(b)
    ss_a = sum(((x-mean_a)**2 for x in a),z)
    ss_b = sum(((x-mean_b)**2 for x in b),z)
    return (ss_a + ss_b) / (len(a) + len(b) - 2)

def dimension_pass(a, b, v, k, counts=None):
    if v == 0:
        return counts[0] == counts[1] if counts is not None else a == b
    return (a-b)**2 <= k*k*v

def pair_pass(a, b, variances, k):
    return all(dimension_pass(a[j],b[j],variances[j],k if j<3 else decimal_fraction(k), (a[4],b[4]) if j==3 else None) for j in range(4))

def exclusion(missing, opposite_raw_n, degree):
    if missing: return REASONS[0]
    if opposite_raw_n == 0: return REASONS[1]
    if degree == 0: return REASONS[2]
    return None

def synthetic_tests():
    out = []
    def add(name, actual, expected):
        assert actual == expected, (name,actual,expected)
        out.append({'test':name,'actual':actual,'expected':expected,'passed':True})
    F=Fraction
    add('zero_variance_equal',dimension_pass(F(2),F(2),F(0),F(1)),True)
    add('zero_variance_unequal',dimension_pass(F(2),F(3),F(0),F(1)),False)
    add('zero_log_variance_uses_original_count_equality',dimension_pass(Decimal(2),Decimal(2),Decimal(0),Decimal(1),(10,11)),False)
    add('exact_caliper_boundary_included',dimension_pass(F(0),F(1),F(1),F(1)),True)
    add('just_outside_boundary_excluded',dimension_pass(F(0),F(1001,1000),F(1),F(1)),False)
    add('half_caliper_boundary_included',dimension_pass(F(0),F(1,2),F(1),F(1,2)),True)
    add('within_group_pooled_ignores_between_group_shift',exact(pooled_variance([F(0),F(2)],[F(10),F(12)])),'2/1')
    add('small_n_scale_is_undefined',pooled_variance([F(0)],[F(1),F(2)]),None)
    a=(F(0),F(0),F(0),Decimal(0),1)
    b1=(F(0),F(2),F(0),Decimal(0),1)
    b2=(F(2),F(0),F(0),Decimal(0),1)
    v=[F(1),F(1),F(1),Decimal(1)]
    add('separate_neighbors_for_separate_dimensions_not_joint',any(pair_pass(a,b,v,F(1)) for b in [b1,b2]),False)
    add('missing_priority_over_no_opposite',exclusion(True,0,0),'missing_invalid_input')
    add('no_opposite_priority_over_no_joint',exclusion(False,0,0),'no_opposite_group_in_prefecture')
    add('no_joint_reason',exclusion(False,1,0),'no_joint_neighbor_within_fixed_caliper')
    add('retained_has_no_exclusion_reason',exclusion(False,1,1),None)
    add('invalid_noninteger_count',allowed_count('1.5'),None)
    add('invalid_negative_count',allowed_count('-1'),None)
    add('blank_is_missing',allowed_count(''),None)
    points=[(F(i,4),F(0),F(0),Decimal(0),1) for i in range(-8,9)]
    s1={i for i,b in enumerate(points) if pair_pass(a,b,v,F(1))}
    s05={i for i,b in enumerate(points) if pair_pass(a,b,v,F(1,2))}
    add('half_support_is_subset',s05 <= s1,True)
    return out
