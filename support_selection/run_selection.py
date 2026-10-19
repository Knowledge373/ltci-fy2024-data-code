#!/usr/bin/env python3
"""Compute fixed support sets from final covariates and monthly classifications."""
import argparse
from collections import Counter
import csv
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import runpy

ROOT = Path(__file__).resolve().parent

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read(path):
    with path.open(encoding='utf-8-sig', newline='') as stream:
        return list(csv.DictReader(stream))

def write(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def ids_sha(ids):
    return hashlib.sha256(json.dumps(sorted(ids), ensure_ascii=False, separators=(',', ':')).encode()).hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    hashes = json.loads((ROOT / 'inputs/input_hashes.json').read_text())
    for relative, expected in hashes.items():
        assert sha(ROOT / relative) == expected, relative
    common = ROOT / 'inputs/common_input_1573.csv'
    institution = ROOT / 'inputs/INSTITUTIONAL_JOIN_1573.csv'
    inst = read(institution)
    assert len(inst) == 1573
    wanted = ['insref_id', 'fiscal_year', 'prefecture_code', 'annual_installments',
              'provisional_none', 'income_low3_share', 'age_75_84_share',
              'age_85plus_share', 'ln_insured_end_fy_count', 'insured_count_end_of_fiscal_year']
    rows = [{k: r[k] for k in wanted} for r in read(common)]
    by = {r['insref_id']: r for r in rows}
    assert len(by) == 1573
    values = {r['insref_id']: (*[Fraction(r[k]) for k in ['income_low3_share', 'age_75_84_share', 'age_85plus_share']],
                             Decimal(r['ln_insured_end_fy_count']), int(r['insured_count_end_of_fiscal_year'])) for r in rows}
    assert all(0 <= v <= 1 for row in values.values() for v in row[:3])
    assert all(x[4] > 0 for x in values.values())
    cids = ['C1_six', 'C2_july_8vs9', 'C4_yes12_vs_yes6']
    groups = {c: {'A': [], 'B': []} for c in cids}
    pref = {}
    months = ['2024-04', '2024-05', '2024-06', '2024-07', '2024-08', '2024-09',
              '2024-10', '2024-11', '2024-12', '2025-01', '2025-02', '2025-03']
    assignments = []
    for r in inst:
        i = r['insref_id']
        assert all(r[k] == by[i][k] for k in ['fiscal_year', 'prefecture_code', 'annual_installments', 'provisional_none'])
        annual = int(r['annual_installments'])
        none = int(r['provisional_none'])
        v = [int(r[k]) for k in months]
        assert sum(v) == annual
        pref[i] = r['prefecture_code']
        assigned = {
            'C1_six': ('A' if none else 'B') if annual == 6 else '',
            'C2_july_8vs9': 'A' if none and v == [0,0,0,1,1,1,1,1,1,1,1,0] else 'B' if none and v == [0,0,0,1,1,1,1,1,1,1,1,1] else '',
            'C4_yes12_vs_yes6': ('A' if annual == 12 else 'B') if not none and annual in {6,12} else ''
        }
        for c, group in assigned.items():
            assignments.append({'contrast_id': c, 'ID': i, 'group': group})
            if group:
                groups[c][group].append(i)
    for gg in groups.values():
        for group in gg:
            gg[group].sort()
    mods = {'early': runpy.run_path(str(ROOT / 'math_early.py')),
            'final_p3': runpy.run_path(str(ROOT / 'math_p3.py'))}
    pins = json.loads((ROOT / 'inputs/set_hashes.json').read_text())['sets']
    results, members, pairs, edges, scales, tests = [], [], [], [], [], []
    with localcontext() as ctx:
        ctx.prec = 80
        for label, mod in mods.items():
            tests.extend(dict(module=label, **row) for row in mod['synthetic_tests']())
        for cid in cids:
            m = mods['final_p3' if cid.startswith('C4') else 'early']
            group = groups[cid]
            variances = [m['pooled_variance']([values[i][j] for i in group['A']], [values[i][j] for i in group['B']]) for j in range(4)]
            for j, variance in enumerate(variances):
                actual = Decimal(variance.numerator) / Decimal(variance.denominator) if j < 3 else variance
                scales.append({'contrast_id': cid, 'variable': m['VARIABLES'][j], 'variance_decimal80': str(actual)})
            for kt, k in m['CALIPERS']:
                pin = next(x for x in pins if x['source_contrast_id'] == cid and x['caliper'] == kt)
                degree = Counter()
                current = set()
                compared = 0
                for aid in group['A']:
                    for bid in group['B']:
                        if pref[aid] != pref[bid]:
                            continue
                        aa, bb = values[aid], values[bid]
                        dims = [m['dimension_pass'](aa[j], bb[j], variances[j], k if j < 3 else m['decimal_fraction'](k), (aa[4], bb[4]) if j == 3 else None) for j in range(4)]
                        accepted = m['pair_pass'](aa, bb, variances, k)
                        compared += 1
                        pairs.append({'analysis': pin['analysis_contrast_id'], 'contrast_id': cid, 'caliper': kt,
                                      'A_ID': aid, 'B_ID': bid, 'pass_dimensions': json.dumps(dims), 'edge': accepted})
                        if accepted:
                            current.add((aid, bid))
                            degree[aid] += 1
                            degree[bid] += 1
                            edges.append({'contrast_id': cid, 'caliper': kt, 'A_ID': aid, 'B_ID': bid})
                retained = {}
                for group_name in ['A', 'B']:
                    opposite = 'B' if group_name == 'A' else 'A'
                    retained[group_name] = []
                    for i in group[group_name]:
                        reason = m['exclusion'](False, sum(pref[j] == pref[i] for j in group[opposite]), degree[i])
                        kept = degree[i] > 0
                        members.append({'analysis': pin['analysis_contrast_id'], 'contrast_id': cid, 'caliper': kt,
                                        'group': group_name, 'ID': i, 'retained': kept, 'neighbor_count': degree[i],
                                        'exclusion_reason': reason})
                        if kept:
                            retained[group_name].append(i)
                all_ids = retained['A'] + retained['B']
                tsha = ids_sha(all_ids)
                assert ids_sha(retained['A']) == pin['source_group_A_IDs_sha256']
                assert ids_sha(retained['B']) == pin['source_group_B_IDs_sha256']
                assert tsha == pin['T_IDs_sha256']
                assert ids_sha(retained[pin['D1_source_group']]) == pin['D1_IDs_sha256']
                assert ids_sha(retained[pin['D0_source_group']]) == pin['D0_IDs_sha256']
                results.append({'analysis': pin['analysis_contrast_id'], 'source_contrast': cid, 'caliper': kt,
                                'n_A': len(retained['A']), 'n_B': len(retained['B']), 'n_T': len(all_ids),
                                'n_D1': len(retained[pin['D1_source_group']]), 'n_D0': len(retained[pin['D0_source_group']]),
                                'prefectures': len({pref[i] for i in all_ids}), 'edges': len(current),
                                'pairs_checked': compared, 'T_IDs_sha256': tsha})
    for filename, data in [('GROUP_ASSIGNMENTS.csv', assignments), ('RETAINED_SETS.csv', results),
                           ('MEMBERSHIP.csv', members), ('PAIR_DECISIONS.csv', pairs),
                           ('EDGES.csv', edges), ('FIXED_SCALES.csv', scales), ('SYNTHETIC_TESTS.csv', tests)]:
        write(out / filename, data)
    print(json.dumps({'status': 'PASS', 'output': str(out), 'precision': 80,
                      'sets': len(results), 'membership_rows': len(members),
                      'same_prefecture_pair_decisions': len(pairs), 'edges': len(edges)}, indent=2))

if __name__ == '__main__':
    main()
