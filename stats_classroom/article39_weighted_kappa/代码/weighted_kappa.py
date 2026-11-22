# 统计课堂39：轻度判成中度，与轻度判成重度，能算同样的不一致吗？
# 数据为教学假设，100个独立患者，每人A/B各一次；同一行两值必须保留为一对。
# 指标为过去24小时疼痛对日常活动的影响，A/B盲评同次就诊资料。
# 无=没有疼痛；轻=有疼痛但不影响活动；中=影响活动但仍能完成；重=部分活动无法完成。
# 在Spyder把项目根设为工作目录，按 # %% 顺序运行；F5可运行全文件。
# UTF-8；依赖numpy、pandas、scikit-learn。结果先在Console查看；可保存到运行结果/Python/。
# draw_article.py用于独立重建文章配图，不是分析脚本必需的绘图步骤。

# %% 1. 读入数据：case区分主例/对照，A/B为0无1轻2中3重，count为人数。
from pathlib import Path
import csv
import os
import sys
import numpy as np
import pandas as pd
from sklearn.metrics import cohen_kappa_score

# 兼容命令行__file__与Spyder分节时没有__file__的环境。
if "__file__" in globals():
    project = Path(__file__).resolve().parents[1]
else:
    project = Path.cwd()
    if project.name == "代码":
        project = project.parent
if not (project / "数据" / "paired_tables.csv").is_file():
    raise FileNotFoundError("请把工作目录设为本篇项目根，再从第1节运行。")
sys.path.insert(0, str(project / "代码"))
from kappa_parts import kappa_parts
with (project / "数据" / "paired_tables.csv").open(encoding="utf-8-sig") as f:
    records = list(csv.DictReader(f))
print("字段示例：", records[0])
case_names = ["main", "control"]
expected_n = 100  # 本题每张表应有100人；统计计算仍使用实际计数之和。
observed_cases = set()
for record in records:
    observed_cases.add(record["case"])
if observed_cases != set(case_names):
    raise ValueError("case只能包含main与control，请检查是否读错文件。")
tables = {}
for case in case_names:
    selected = [r for r in records if r["case"] == case]  # r['case']取字段，条件筛选主例或对照。
    if len(selected) != 16:
        raise ValueError("每张方表须有16格，包括零格。")
    tab = np.zeros((4, 4), dtype=int)
    seen = set()
    for r in selected:
        i = int(r["A"])
        j = int(r["B"])
        count = int(r["count"])
        if not (0 <= i <= 3 and 0 <= j <= 3 and count >= 0) or (i, j) in seen:
            raise ValueError("等级、计数或重复格有误。")
        seen.add((i, j))
        tab[i, j] = count  # NumPy位置从0开始，正好对应CSV的A、B编码。
    if tab.sum() != expected_n:
        raise ValueError(f"{case}总人数应为{expected_n}，实际为{tab.sum()}。")
    tables[case] = tab
    print(case, "行A列B：无/轻/中/重\n", tab)

# 日常练习先只看Console；在此处把这个布尔值改为True即可保存结果。
save_results = False

# %% 2. 用同一顺序构造一致权重；所有格按原位置相乘，不改行列顺序。
levels = np.arange(4)
# [:, None]把等级排成一列，再与一行等级相减，得到16格的等级差。
distance = np.abs(levels[:, None] - levels[None, :])
weights = {}
# eye(4)的对角为1、其余为0，对应未加权的完全同级记分。
weights["unweighted"] = np.eye(4)
weights["linear"] = 1 - distance / 3
weights["quadratic"] = 1 - (distance / 3) ** 2
print("线性一致权重：\n", weights["linear"])
print("二次一致权重：\n", weights["quadratic"])

# %% 3.1 主例线性κ：把正文手算逐步写成代码。
main_manual = tables["main"]
manual_n = main_manual.sum()
manual_exact = np.trace(main_manual) / manual_n  # trace取主对角线，即A/B完全同级格。
manual_p = main_manual / manual_n
manual_row = manual_p.sum(axis=1)  # axis=1横向求和，得到A行边际。
manual_col = manual_p.sum(axis=0)  # axis=0纵向求和，得到B列边际。
manual_chance = np.outer(manual_row, manual_col)  # 外积构造A、B独立时的16格参照比例。
manual_weight = weights["linear"]
manual_po = np.sum(manual_weight * manual_p)  # 逐格相乘一致权重与真实比例，再求和。
manual_pe = np.sum(manual_weight * manual_chance)  # 同样逐格加权独立参照，再求和。
manual_kappa = (manual_po - manual_pe) / (1 - manual_pe)
distance_steps = []
for d in range(4):
    mask = distance == d  # mask筛选同距离格；同距离格共享同一权重。
    patients = int(main_manual[mask].sum())
    chance_cells = int(mask.sum())
    chance_probability = float(manual_chance[mask].sum())
    linear_weight = float(manual_weight[mask][0])
    Po_term = linear_weight * patients / manual_n  # patients是该距离的患者人数。
    Pe_term = linear_weight * chance_probability  # chance_probability是边际外积参照比例之和。
    distance_steps.append([d, patients, chance_cells, chance_probability, linear_weight, Po_term, Pe_term])
distance_steps = pd.DataFrame(distance_steps, columns=["distance", "patients", "chance_cells",
    "chance_probability", "linear_weight", "Po_term", "Pe_term"])
assert np.isclose(distance_steps["Po_term"].sum(), manual_po, atol=1e-12, rtol=0)
assert np.isclose(distance_steps["Pe_term"].sum(), manual_pe, atol=1e-12, rtol=0)
linear_steps = pd.DataFrame([{
    "n": manual_n,
    "exact_agreement": manual_exact,
    "Po": manual_po,
    "Pe": manual_pe,
    "kappa": manual_kappa,
}])
print(f"3.1主例：n={manual_n}；完全同级比例={manual_exact:.2%}（回答有多少人两次评级相同）。")
print("每格联合比例p：\n", manual_p)
print("行边际r与列边际c：\n", manual_row, "\n", manual_col)
print("独立机会参照r×c：\n", manual_chance)
print(f"线性Po={manual_po:.7f}（真实配对的加权一致量）；Pe={manual_pe:.7f}（边际独立时的参照量）；κ={manual_kappa:.7f}（相对机会基准的加权一致性）。")
print("linear_steps（3.1逐步结果）：\n", linear_steps)
print("distance_steps（按距离分解式4/5）：\n", distance_steps.to_string(index=False))

# %% 3.2 全表比较、一般SE/Wald及原生核对。
# kappa_parts是随项目提供的自编函数，返回Po、Pe、κ、SE、Wald下/上限。
control = tables["control"]
distance_comparison_rows = []
for d in range(4):
    mask = distance == d
    distance_comparison_rows.append([d, int(main_manual[mask].sum()), int(control[mask].sum())])
distance_comparison = pd.DataFrame(distance_comparison_rows, columns=["distance", "main_count", "control_count"])
control_checks = pd.DataFrame([{
    "same_diagonal": np.array_equal(np.diag(main_manual), np.diag(control)),
    "same_row_totals": np.array_equal(main_manual.sum(axis=1), control.sum(axis=1)),
    "same_column_totals": np.array_equal(main_manual.sum(axis=0), control.sum(axis=0)),
}])
print("distance_comparison（主例/对照按距离计数）：\n", distance_comparison)
print("control_checks（对照保持不变条件）：\n", control_checks)
if not bool(control_checks.iloc[0].all()):
    raise ValueError("control应保持主例的对角和行列合计，请检查分歧移动。")

result_rows = []
for case, tab in tables.items():
    exact_agreement = np.trace(tab) / tab.sum()
    print(f"{case}完全同级比例：{exact_agreement:.2%}（对角人数/总人数）。")
    # 按每格人数展开：相同位置的A/B是同一患者的两次评级。
    paired = []
    for i in range(4):
        for j in range(4):
            paired.extend([(i, j)] * int(tab[i, j]))  # tuple保留同一患者的A/B联合位置，并按人数重复。
    paired = np.asarray(paired, dtype=int)
    # 这个循环依次覆盖2张表×3种权重；每次都从当前表计算，不使用手抄结果。
    for method, w in weights.items():
        values = kappa_parts(tab, w)  # values[0..5]依次是Po、Pe、κ、SE、Wald下限、上限。
        # values按Po、Pe、κ、SE、CI下限、CI上限排列，Python位置从0开始。
        observed_loss = 1 - values[0]
        expected_loss = 1 - values[1]
        loss_ratio = np.nan
        if expected_loss > 1e-14:
            loss_ratio = observed_loss / expected_loss
        # None表示未加权；labels固定全部4个等级，即使某级本次没有出现。
        native_weight = None if method == "unweighted" else method
        native = np.nan
        if np.isfinite(values[2]):
            # [:, 0]取全部患者的A评级，[:, 1]取同一顺序的B评级。
            native = cohen_kappa_score(paired[:, 0], paired[:, 1], labels=[0, 1, 2, 3], weights=native_weight)
            np.testing.assert_allclose(values[2], native, atol=1e-12, rtol=0)
        result_rows.append([case, method, *values, exact_agreement, observed_loss, expected_loss, loss_ratio, native])
        print(case, method, "Po/Pe/κ/SE/CI：", np.round(values, 6), "sklearn κ=", round(native, 6))
        print(f"实际平均损失Lo={observed_loss:.6f}，机会平均损失Le={expected_loss:.6f}，Lo/Le={loss_ratio:.6f}。")
print("κ按指定规则比较一致程度；完全同级比例另由对角人数除以总人数计算。")
result_columns = ["case", "method", "Po", "Pe", "kappa", "SE", "CI_low", "CI_high",
                  "exact_agreement", "Lo", "Le", "loss_ratio", "sklearn_kappa"]
results = pd.DataFrame(result_rows, columns=result_columns)
# &把两个条件连接起来；iloc[0]取筛选后的首行，得到主例线性结果。
main_linear = results[(results["case"] == "main") & (results["method"] == "linear")].iloc[0]
main_quadratic = results[(results["case"] == "main") & (results["method"] == "quadratic")].iloc[0]
np.testing.assert_allclose(linear_steps.loc[0, ["Po", "Pe", "kappa"]].to_numpy(dtype=float),
                           main_linear[["Po", "Pe", "kappa"]].to_numpy(dtype=float), atol=1e-12, rtol=0)
function_comparison = results[["case", "method", "kappa", "sklearn_kappa"]].rename(
    columns={"kappa": "kappa_formula", "sklearn_kappa": "kappa_native"})
print("function_comparison（自编公式与sklearn原生点估计）：\n", function_comparison)
# report_text为由当前结果格式化生成的报告句，不调用新的report软件。
report_text = (f"对{manual_n}名患者按无、轻、中、重顺序分析；预设线性一致权重为1−|i−j|/3。"
               f"完全同级比例为{main_linear['exact_agreement']:.1%}，线性Po={main_linear['Po']:.4f}、Pe={main_linear['Pe']:.4f}、"
               f"κ={main_linear['kappa']:.3f}，95%Wald置信区间为[{main_linear['CI_low']:.3f}, {main_linear['CI_high']:.3f}]；"
               f"二次敏感性分析κ={main_quadratic['kappa']:.3f}，95%Wald置信区间为"
               f"[{main_quadratic['CI_low']:.3f}, {main_quadratic['CI_high']:.3f}]。")
print("report_text：", report_text)

# %% 4. 成对Bootstrap：固定权重，每轮重抽100条完整患者配对记录。
# 只分析主例；两种权重使用同一批重抽表，避免额外引入不同的模拟样本。
tab = tables["main"]
pair_cells = np.repeat(np.arange(16), tab.ravel())
# pair_cells中每项代表一个人的联合格编号；ravel按行展开，bincount后同样按行还原。
repetitions = 10000
rng = np.random.default_rng(20261122)
boot = np.full((repetitions, 2), np.nan)
for b in range(repetitions):
    # size等于患者数，replace=True表示有放回；抽的是完整联合格记录。
    chosen = rng.choice(pair_cells, size=len(pair_cells), replace=True)  # replace=True表示有放回重抽。
    sampled = np.bincount(chosen, minlength=16).reshape(4, 4)  # bincount计数后按同一行展开顺序重构4×4表。
    boot[b, 0] = kappa_parts(sampled, weights["linear"])[2]
    boot[b, 1] = kappa_parts(sampled, weights["quadratic"])[2]
bootstrap_rows = []
for j, method in enumerate(["linear", "quadratic"]):
    valid = np.isfinite(boot[:, j])
    if not np.any(valid):
        raise ValueError("全部重抽κ无定义，请检查数据与等级。")
    # 0.025、0.975为双侧95%百分位区间的两个端点概率。
    # method=linear与R type=7同为位置1+(B-1)q上的线性插值。
    interval = np.quantile(boot[valid, j], [0.025, 0.975], method="linear")
    row = [method, repetitions, int(valid.sum()), int((~valid).sum()), *interval]
    bootstrap_rows.append(row)
    print(f"主例{method}：{repetitions}轮中有效{valid.sum()}轮，无定义{(~valid).sum()}轮；95%百分位区间[{interval[0]:.6f}, {interval[1]:.6f}]。")
bootstrap_columns = ["method", "B", "valid", "undefined", "CI_low", "CI_high"]
bootstrap_results = pd.DataFrame(bootstrap_rows, columns=bootstrap_columns)
print("bootstrap_results（可在Variable Explorer查看）：\n", bootstrap_results)
if any(row[3] > 0 for row in bootstrap_rows):
    print("有重抽轮次的κ无定义，请结合报告的次数检查这份数据的重抽分布。")

# %% 5. 可选保存；首次练习先在Console查看，环境变量只用于批量复现输出。
# 日常只改第1节的save_results；环境变量1用于批量复现时覆盖为True。
if os.environ.get("ARTICLE39_SAVE_RESULTS") == "1":
    save_results = True
if save_results:
    out = project / "运行结果" / "Python"
    out.mkdir(parents=True, exist_ok=True)
    with (out / "summary.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["case", "method", "Po", "Pe", "kappa", "SE", "CI_low", "CI_high", "exact_agreement", "Lo", "Le", "loss_ratio", "sklearn_kappa"])
        writer.writerows(result_rows)
    with (out / "bootstrap.csv").open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["method", "B", "valid", "undefined", "CI_low", "CI_high"])
        writer.writerows(bootstrap_rows)
# 练习：先预测把近处分歧移远的影响，再读control结果；不按κ最高来挑权重。
