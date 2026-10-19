#!/usr/bin/env python3
"""Run the fixed exploratory calculations in a new output directory."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    hashes = json.loads((ROOT / 'inputs/input_hashes.json').read_text())
    for relative, expected in hashes.items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected, relative
    inp = ROOT / 'inputs/common_input_1573.csv'
    base = ROOT / 'inputs/baseline'
    digest = hashlib.sha256(inp.read_bytes()).hexdigest()
    assert digest == 'ea25bf07ed1ff42626ca73719ca252af7f563ca35cd19029139863cdbca468f2'
    for name in ['weights', 'influence_distribution', 'institution_size', 'logs']:
        (out / name).mkdir()
    cfg = {
        'input_path': str(inp), 'input_sha256': digest,
        'baseline_path': str(base / 'analysis_rows_1573.csv'),
        'amount_units': 'JPY, FY2024 current-assessment flow',
        'columns': {
            'N': 'insured_count_end_of_fiscal_year',
            'A': 'ordinary_current_assessment_yen',
            'C': 'ordinary_current_collected_yen',
            'S': 'special_current_assessment_yen',
            'T': 'total_current_assessment_yen'
        }
    }
    config = out / 'institution_size/config.json'
    config.write_text(json.dumps(cfg, indent=2) + '\n')
    env = os.environ.copy()
    env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', MPLCONFIGDIR=str(out / 'matplotlib_cache'),
               PYTHONDONTWRITEBYTECODE='1')
    jobs = [
        ('weights', 'weights/analyze_weights.py', [
            '--input', inp, '--expected-sha256', digest,
            '--count-column', 'insured_count_end_of_fiscal_year',
            '--ordinary-assessment-column', 'ordinary_current_assessment_yen',
            '--total-assessment-column', 'total_current_assessment_yen',
            '--ordinary-collection-column', 'ordinary_current_collected_yen']),
        ('weights', 'weights/annotate_contrast_support.py', []),
        ('influence_distribution', 'influence_distribution/analyze_distribution_influence.py', [
            '--baseline', base, '--extended', inp, '--extended-sha256', digest,
            '--assessed-column', 'ordinary_current_assessment_yen',
            '--received-column', 'ordinary_current_collected_yen',
            '--population-column', 'insured_count_end_of_fiscal_year', '--amount-unit', 'yen',
            '--out', out / 'influence_distribution']),
        ('institution_size', 'institution_size/analyze.py', ['--config', config]),
        ('external_work', 'external_work/explore.py', [
            '--input', inp, '--model-specification', base / 'model_specification.json',
            '--out', out / 'external_work'])
    ]
    for name, script, argv in jobs:
        child_env = env.copy()
        child_env['EXPLORATION_OUTPUT_DIR'] = str(out / name)
        with (out / 'logs' / f'{name}_{Path(script).stem}.stdout.txt').open('w') as stdout, \
             (out / 'logs' / f'{name}_{Path(script).stem}.stderr.txt').open('w') as stderr:
            subprocess.run([sys.executable, str(ROOT / script), *map(str, argv)],
                           env=child_env, stdout=stdout, stderr=stderr, check=True)
        print(f'{name}: PASS', flush=True)
    print(json.dumps({'status': 'PASS', 'output': str(out)}))

if __name__ == '__main__':
    main()
