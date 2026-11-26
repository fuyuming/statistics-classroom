"""与car::qqPlot(fit, distribution='t', line='quartiles', envelope=FALSE)对应。"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats


def studentized_qq(fit, draw=True):
    z = np.sort(np.asarray(fit.get_influence().resid_studentized_external))
    n = len(z)
    degrees = float(fit.df_resid-1)
    if degrees <= 0 or not np.all(np.isfinite(z)):
        raise ValueError("finite externally studentized residuals and df_resid > 1 required")
    a = 0.375 if n <= 10 else 0.5
    probability = (np.arange(1,n+1)-a)/(n+1-2*a)
    theoretical = stats.t.ppf(probability, df=degrees)
    sample_quartiles = np.quantile(z,[0.25,0.75],method="linear")
    theory_quartiles = stats.t.ppf([0.25,0.75],df=degrees)
    slope = np.diff(sample_quartiles)[0]/np.diff(theory_quartiles)[0]
    intercept = sample_quartiles[0]-slope*theory_quartiles[0]
    points = pd.DataFrame({"probability":probability,"theoretical":theoretical,
                           "studentized":z,"reference_line":intercept+slope*theoretical})
    out = {"points":points,"df":degrees,"slope":slope,"intercept":intercept}
    if draw:
        figure, axes = plt.subplots(figsize=(7,5))
        axes.scatter(theoretical,z,color="0.2")
        axes.plot(theoretical,points["reference_line"],color="blue")
        axes.set(xlabel=f"Theoretical t quantiles (df={degrees:g})",
                 ylabel="Externally studentized residual",title="Studentized residual Q-Q plot")
        figure.tight_layout()
        out.update(figure=figure,axes=axes)
    return out
