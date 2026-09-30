"""跨语言核对：需先运行Python、R、MATLAB主脚本。不是独立性统计检验。"""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parents[1]
# 公开运行日志仅清理显示用的行尾空格，不改变任何数值或内容。
for log in (root / '运行结果').glob('*/运行记录.txt'):
    lines = log.read_text(encoding='utf-8').splitlines()
    log.write_text('\n'.join(line.rstrip() for line in lines) + '\n', encoding='utf-8')
files = ['position_results.csv', 'covariance_identity.csv', 'sampling_summary.csv',
         'conditional_summary.csv', 'pair_summary.csv', 'figure2_coordinates.csv',
         'normal_statistics.csv', 'shifted_exponential_statistics.csv']
rows = []
for name in files:
    baseline = pd.read_csv(root / '运行结果' / 'Python' / name)
    for language in ['R', 'MATLAB']:
        candidate = pd.read_csv(root / '运行结果' / language / name)
        assert list(candidate.columns) == list(baseline.columns), (name, language, 'columns')
        assert candidate.shape == baseline.shape
        numeric = baseline.select_dtypes(include='number').columns
        for column in baseline.columns.difference(numeric):
            assert baseline[column].equals(candidate[column])
        left = baseline[numeric].to_numpy()
        right = candidate[numeric].to_numpy()
        # 同一输入应仅有浮点和CSV保存精度差；不允许用模拟随机误差掩盖计算差异。
        assert np.allclose(left, right, atol=1e-10, rtol=1e-10), (name, language)
        rows.append({'file': name, 'compared_with_python': language,
                     'rows': len(baseline), 'max_absolute_difference': np.max(np.abs(left-right)),
                     'passed': True})
checks = pd.DataFrame(rows)
checks.to_csv(root / '运行结果' / 'cross_language_checks.csv', index=False)
hash_lines = []
for path in sorted((root / '数据').glob('*')):
    hash_lines.append(hashlib.sha256(path.read_bytes()).hexdigest() + '  ' + path.name)
(root / '数据' / 'SHA256SUMS.txt').write_text(
    '\n'.join(line for line in hash_lines if not line.endswith('SHA256SUMS.txt')) + '\n', encoding='utf-8')
print(checks.to_string(index=False))
