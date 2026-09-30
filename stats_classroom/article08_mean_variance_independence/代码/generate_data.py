"""生成本课共享模拟输入。仅重新抽样时运行，普通练习直接读取已有 CSV。
数据完全由指定模型产生，不是真实实验数据。UTF-8，无缺失值，无量纲。
种子20260930；NumPy default_rng/PCG64。三语言读取同一文件以便逐项核对。
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
folder = root / '数据'
folder.mkdir(exist_ok=True)
rng = np.random.default_rng(20260930)
repetitions = 20000
sample_size = 10

# 每行是一份独立样本，x1至x10为这份样本中的10个独立观测。
# 两种总体均值都为10、方差都为4，差别在分布形状。
for model in ['normal', 'shifted_exponential']:
    if model == 'normal':
        values = rng.normal(loc=10, scale=2, size=(repetitions, sample_size))
    else:
        values = 8 + rng.exponential(scale=2, size=(repetitions, sample_size))
    frame = pd.DataFrame(values, columns=[f'x{i}' for i in range(1, 11)])
    frame.insert(0, 'sample_id', np.arange(1, repetitions + 1))
    frame.to_csv(folder / f'{model}_samples.csv', index=False, float_format='%.17g')

# 图2的三种模型共享X；independent_y独立生成；sign独立取正负1。
x = rng.normal(size=1200)
independent_y = rng.normal(size=1200)
sign = rng.choice([-1, 1], size=1200)
pairs = pd.DataFrame({'pair_id': np.arange(1, 1201), 'x': x,
                      'independent_y': independent_y, 'sign': sign})
pairs.to_csv(folder / 'joint_pairs.csv', index=False, float_format='%.17g')
examples = pd.DataFrame({'group': ['A', 'B', 'C'], 'x1': [8, 18, 6],
                         'x2': [10, 20, 10], 'x3': [12, 22, 14]})
examples.to_csv(folder / 'position_examples.csv', index=False)
metadata = {'seed': 20260930, 'generator': 'NumPy default_rng / PCG64',
            'numpy': np.__version__, 'repetitions': repetitions, 'n': sample_size,
            'mu': 10, 'sigma_squared': 4, 'normal': 'N(10,4)',
            'shifted_exponential': '8 + Exponential(scale=2)',
            'pair_count': 1200, 'missing_values': 0,
            'units': 'dimensionless', 'source': 'synthetic teaching data'}
(folder / 'generation_parameters.json').write_text(
    json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
print('已生成共享教学数据。重新生成后应重跑三套计算及配图。')
