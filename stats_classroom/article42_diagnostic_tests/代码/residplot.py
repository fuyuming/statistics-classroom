"""按R in Action残差分布图思路编写的跨语言课堂函数。"""
import numpy as np
import matplotlib.pyplot as plt
from statsmodels.regression.linear_model import OLS


def residplot(fit, nbreaks=10, draw=True):
    """输入普通OLS拟合结果；返回外部学生化残差、直方图、曲线及带宽。

    nbreaks是等宽区间数。draw=False时只计算，不创建图。
    """
    if not isinstance(nbreaks, (int, np.integer)) or isinstance(nbreaks, bool) or nbreaks < 2:
        raise ValueError("nbreaks must be an integer >= 2")
    if not hasattr(fit, "model") or not isinstance(fit.model, OLS):
        raise ValueError("fit must be a statsmodels OLS result")
    design = np.asarray(fit.model.exog)
    if np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("the design matrix must have full column rank")
    # 外部学生化采用删除当前观测后的误差方差；不要换成内部标准化残差。
    z = np.asarray(fit.get_influence().resid_studentized_external, dtype=float)
    n = len(z)
    if fit.df_resid <= 1 or not np.all(np.isfinite(z)) or np.ptp(z) == 0:
        raise ValueError("finite, varying studentized residuals and df_resid > 1 are required")
    sample_sd = np.std(z, ddof=1)
    q1, q3 = np.quantile(z, [0.25, 0.75], method="linear")
    spread = min(sample_sd, (q3-q1)/1.34)
    if spread <= 0:
        spread = sample_sd
    bandwidth = 0.9*spread*n**(-0.2)
    edges = np.linspace(z.min(), z.max(), nbreaks+1)
    counts, _ = np.histogram(z, bins=edges)
    widths = np.diff(edges)
    density = counts/(n*widths)
    grid = np.linspace(z.min()-3*bandwidth, z.max()+3*bandwidth, 512)
    normal = np.exp(-0.5*((grid-z.mean())/sample_sd)**2)/(sample_sd*np.sqrt(2*np.pi))
    # 显式共用Gaussian核和带宽，避免各语言默认平滑参数不同。
    distances = (grid[:, None]-z[None, :])/bandwidth
    kde = np.mean(np.exp(-0.5*distances**2), axis=1)/(bandwidth*np.sqrt(2*np.pi))
    out = {"studentized_residuals": z,
           "histogram": {"left": edges[:-1], "right": edges[1:], "count": counts, "density": density},
           "curves": {"x": grid, "normal": normal, "kde": kde}, "bandwidth": float(bandwidth)}
    if draw:
        figure, axes = plt.subplots(figsize=(8, 5))
        axes.bar(edges[:-1], density, width=widths, align="edge", color="0.78", edgecolor="white")
        axes.plot(grid, normal, color="blue", label="Fitted normal")
        axes.plot(grid, kde, "r--", label="Gaussian KDE")
        axes.plot(z, np.zeros_like(z), "|", color="brown", markersize=12, label="Residual positions")
        axes.set_xlim(grid[0], grid[-1])
        axes.set_ylim(0, 1.08*max(density.max(), normal.max(), kde.max()))
        axes.set(xlabel="Externally studentized residual", ylabel="Density", title="Residual distribution")
        axes.legend()
        figure.tight_layout()
        out.update(figure=figure, axes=axes)
    return out
