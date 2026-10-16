#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""四个完全假设数据的确定性教学算例。

本脚本只做枚举、加权平均和普通最小二乘（OLS）。它没有调用 quap，
也没有重复模拟，因此输出用于课堂算术演示，不是贝叶斯后验或原课复现。
"""

from __future__ import annotations

import csv
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "input_data.csv"
OUTPUT = ROOT / "results" / "python_results.csv"


def read_rows() -> list[dict[str, str]]:
    with INPUT.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in (None, "") else math.nan


def weighted_mean(rows: list[dict[str, str]], value_key: str, weight_key: str) -> float:
    weights = [number(row, weight_key) for row in rows]
    values = [number(row, value_key) for row in rows]
    return sum(value * weight for value, weight in zip(values, weights)) / sum(weights)


def solve_linear(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """用带主元的高斯消元解小型正规方程，避免引入数值计算依赖。"""
    n = len(vector)
    augmented = [row[:] + [rhs] for row, rhs in zip(matrix, vector)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda row: abs(augmented[row][col]))
        if abs(augmented[pivot][col]) < 1e-14:
            raise ValueError("设计矩阵不可逆")
        augmented[col], augmented[pivot] = augmented[pivot], augmented[col]
        divisor = augmented[col][col]
        augmented[col] = [entry / divisor for entry in augmented[col]]
        for row in range(n):
            if row == col:
                continue
            factor = augmented[row][col]
            augmented[row] = [
                left - factor * right
                for left, right in zip(augmented[row], augmented[col])
            ]
    return [augmented[row][-1] for row in range(n)]


def ols_coefficients(x_columns: list[list[float]], y: list[float]) -> list[float]:
    """返回含截距的 OLS 系数；x_columns 是各自变量列。"""
    design = [[1.0] + [column[row] for column in x_columns] for row in range(len(y))]
    xtx = [
        [sum(design[row][i] * design[row][j] for row in range(len(y))) for j in range(len(design[0]))]
        for i in range(len(design[0]))
    ]
    xty = [sum(design[row][i] * y[row] for row in range(len(y))) for i in range(len(design[0]))]
    return solve_linear(xtx, xty)


def metric(rows: list[tuple[int, str, float]]) -> dict[str, float | int | str]:
    example, name, value = rows[-1]
    return {"example": example, "metric": name, "value": value}


def main() -> None:
    rows = read_rows()
    results: list[dict[str, float | int | str]] = []

    # 算例 1：先按真实分层比例求原始均值，再把高、低分档各设为 50%。
    ex1 = [row for row in rows if row["example"] == "1"]
    trained = [row for row in ex1 if row["trained"] == "1"]
    untrained = [row for row in ex1 if row["trained"] == "0"]
    raw_trained = weighted_mean(trained, "mean_value", "count")
    raw_untrained = weighted_mean(untrained, "mean_value", "count")
    std_trained = sum(number(row, "mean_value") for row in ex1 if row["trained"] == "1") / 2
    std_untrained = sum(number(row, "mean_value") for row in ex1 if row["trained"] == "0") / 2
    results.extend(
        [
            metric([(1, "raw_trained_mean", raw_trained)]),
            metric([(1, "raw_untrained_mean", raw_untrained)]),
            metric([(1, "raw_difference", raw_trained - raw_untrained)]),
            metric([(1, "standardized_trained_mean", std_trained)]),
            metric([(1, "standardized_untrained_mean", std_untrained)]),
            metric([(1, "standardized_difference", std_trained - std_untrained)]),
        ]
    )

    # 算例 2：条件概率 P(乙运动 | 甲运动) 与 P(乙运动 | 甲不运动)。
    ex2 = [row for row in rows if row["example"] == "2"]
    for phase in ("before", "after"):
        cells = {row["friend_cell"]: number(row, "count") for row in ex2 if row["phase"] == phase}
        p_b_given_a = cells["both"] / (cells["both"] + cells["a_only"])
        p_b_given_not_a = cells["b_only"] / (cells["b_only"] + cells["neither"])
        results.extend(
            [
                metric([(2, f"{phase}_P_B_exercises_given_A_exercises", p_b_given_a)]),
                metric([(2, f"{phase}_P_B_exercises_given_A_not_exercises", p_b_given_not_a)]),
            ]
        )

    # 算例 3：Y~X 的斜率为 1/130；加入 Z 后斜率为 0.2。
    # U 是造成 Y 差异的变量；这里只报告有限四行的算术系数，不报告显著性或 SE。
    ex3 = [row for row in rows if row["example"] == "3"]
    x = [number(row, "x") for row in ex3]
    y = [number(row, "y") for row in ex3]
    z = [number(row, "z") for row in ex3]
    coef_unadjusted = ols_coefficients([x], y)
    coef_adjusted = ols_coefficients([x, z], y)
    results.extend(
        [
            metric([(3, "ols_y_on_x_x_coefficient", coef_unadjusted[1])]),
            metric([(3, "ols_y_on_x_plus_z_x_coefficient", coef_adjusted[1])]),
        ]
    )

    # 算例 4：按培训率混合两种培训状态，再比较固定培训状态下的经验差。
    ex4 = [row for row in rows if row["example"] == "4"]
    by_year: dict[int, list[dict[str, str]]] = {}
    for row in ex4:
        by_year.setdefault(int(number(row, "experience_years")), []).append(row)
    overall: dict[int, float] = {}
    for year, year_rows in by_year.items():
        rate = number(year_rows[0], "training_rate")
        weights = [(rate if row["trained"] == "1" else 1 - rate) for row in year_rows]
        overall[year] = sum(number(row, "mean_value") * weight for row, weight in zip(year_rows, weights))
    untrained_by_year = {int(number(row, "experience_years")): number(row, "mean_value") for row in ex4 if row["trained"] == "0"}
    trained_by_year = {int(number(row, "experience_years")): number(row, "mean_value") for row in ex4 if row["trained"] == "1"}
    results.extend(
        [
            metric([(4, "year1_overall_mean", overall[1])]),
            metric([(4, "year2_overall_mean", overall[2])]),
            metric([(4, "overall_difference", overall[2] - overall[1])]),
            metric([(4, "experience_difference_fixed_untrained", untrained_by_year[2] - untrained_by_year[1])]),
            metric([(4, "experience_difference_fixed_trained", trained_by_year[2] - trained_by_year[1])]),
            metric([(4, "training_difference_year1", trained_by_year[1] - untrained_by_year[1])]),
            metric([(4, "training_difference_year2", trained_by_year[2] - untrained_by_year[2])]),
        ]
    )

    expected = {
        "raw_difference": 28,
        "standardized_difference": 10,
        "before_P_B_exercises_given_A_exercises": 0.5,
        "before_P_B_exercises_given_A_not_exercises": 0.5,
        "after_P_B_exercises_given_A_exercises": 0.8,
        "after_P_B_exercises_given_A_not_exercises": 0.2,
        "ols_y_on_x_x_coefficient": 1 / 130,
        "ols_y_on_x_plus_z_x_coefficient": 0.2,
        "overall_difference": 6,
        "experience_difference_fixed_untrained": 2,
        "experience_difference_fixed_trained": 2,
    }
    for result in results:
        if result["metric"] in expected:
            assert math.isclose(float(result["value"]), expected[result["metric"]], rel_tol=1e-10, abs_tol=1e-10)

    OUTPUT.parent.mkdir(exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["example", "metric", "value"])
        writer.writeheader()
        writer.writerows(results)
    print("Python：四个算例的预期数字核对通过。")
    print(f"结果已写入：{OUTPUT}")


if __name__ == "__main__":
    main()
