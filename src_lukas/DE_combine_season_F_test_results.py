import pandas as pd
import numpy as np
import yaml
from misc import de_fn
import sys

config_fn = sys.argv[1]
with open(config_fn, 'r') as f:
    config = yaml.load(f, Loader=yaml.loader.FullLoader)
    model = config['model']
    celltype = config['celltype']
    results_dir = config['results_dir']

steal_suffix = '_V2_V3'
new_coefs = {
    'SEASON.V2': '_DONOR.IC_DATE_2PI_SIN.V2_V2.DONOR.IC_DATE_2PI_COS',
    'SEASON.V3': '_DONOR.IC_DATE_2PI_SIN.V3_V3.DONOR.IC_DATE_2PI_COS',
    'SEASON_correction': '_SAMPLE.VISIT_DATE_2PI_SIN_SAMPLE.VISIT_DATE_2PI_COS'
}

season_df = pd.read_csv(de_fn(celltype, model, results_dir=results_dir).replace(
    '.csv', '{}.csv'.format(steal_suffix if steal_suffix else '')), index_col=0)

for coef in new_coefs:
    df = pd.read_csv(de_fn(celltype, model, results_dir=results_dir).replace(
        '.csv', '{}.csv'.format(new_coefs[coef] if new_coefs[coef] else '')), index_col=0)
    assert df.index.equals(season_df.index)
    season_df['F.{}'.format(coef)] = df['F']
    season_df['p.value.{}'.format(coef)] = df['F.p.value']

season_df.to_csv(de_fn(celltype, model, results_dir=results_dir))

new_contrasts = [list(new_coefs.keys())]
config['contrasts'] = new_contrasts
config['enr_contrasts'] = new_contrasts
with open(config_fn, 'w') as f:
    yaml.dump(config, f, default_flow_style=False, sort_keys=False)

print(model)
print(season_df.head(), '\n')
