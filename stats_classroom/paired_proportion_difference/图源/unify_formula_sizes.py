"""Re-render the five recent formula assets at the reference 20pt scale."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties, findfont

BASE = Path(__file__).resolve().parents[1] / "文章配图"
def find_chinese_font():
    for name in ("Noto Sans CJK SC", "Microsoft YaHei", "SimHei", "PingFang SC", "WenQuanYi Micro Hei"):
        try:
            return findfont(FontProperties(family=name), fallback_to_default=False)
        except ValueError:
            continue
    raise RuntimeError("请先安装可用中文字体，如思源黑体、黑体或微软雅黑。")

FONT = find_chinese_font()
BLUE, INK, MUTED, RED = "#0000F5", "#243746", "#657A7D", "#D9534F"
CN = FontProperties(fname=FONT)
W, DPI = 7.2, 220


def new_canvas(height):
    fig, ax = plt.subplots(figsize=(W, height), dpi=DPI)
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    return fig, ax


def math(ax, expr, y, size=20, x=.5):
    return ax.text(x, y, "$" + expr + "$", transform=ax.transAxes,
                   ha="center", va="center", fontsize=size, color=INK)


def note(ax, text, y, x=.5, color=MUTED):
    ax.text(x, y, text, transform=ax.transAxes, ha="center", va="center",
            fontsize=13, color=color, fontproperties=CN)


def save(fig, name):
    path = BASE / name
    fig.savefig(path, dpi=DPI, facecolor="white")  # deliberately no tight bbox
    plt.close(fig)
    print(path)


def formula03_cancellation():
    fig, ax = new_canvas(5.3)
    ax.text(.05, .96, "式（4）的系数推导", transform=ax.transAxes, ha="left", va="center",
            fontsize=13, color=BLUE, fontproperties=CN, fontweight="bold")
    ax.text(.95, .96, "s = n − r − t", transform=ax.transAxes, ha="right", va="center",
            fontsize=13, color=INK)
    math(ax, r"C(n,r)\,\times\,C(n-r,t)", .81)
    math(ax, r"=\left[\frac{n!}{r!(n-r)!}\right]\times\left[\frac{(n-r)!}{t!s!}\right]", .66)
    ax.text(.10, .49, "=", transform=ax.transAxes, ha="left", va="center", fontsize=20, color=INK)
    nums = [ax.text(.34, .55, r"$n!$", transform=ax.transAxes, ha="center", va="center", fontsize=20, color=INK),
            ax.text(.54, .55, r"$(n-r)!$", transform=ax.transAxes, ha="center", va="center", fontsize=20, color=INK)]
    dens = [ax.text(.27, .43, r"$r!$", transform=ax.transAxes, ha="center", va="center", fontsize=20, color=INK),
            ax.text(.47, .43, r"$(n-r)!$", transform=ax.transAxes, ha="center", va="center", fontsize=20, color=INK),
            ax.text(.67, .43, r"$t!$", transform=ax.transAxes, ha="center", va="center", fontsize=20, color=INK),
            ax.text(.77, .43, r"$s!$", transform=ax.transAxes, ha="center", va="center", fontsize=20, color=INK)]
    ax.plot([.20, .84], [.49, .49], color=INK, lw=1.5, transform=ax.transAxes)
    fig.canvas.draw(); renderer = fig.canvas.get_renderer()
    for txt in (nums[1], dens[1]):
        bb = txt.get_window_extent(renderer=renderer).transformed(ax.transAxes.inverted())
        ax.plot([bb.x0, bb.x1], [bb.y0, bb.y1], color=RED, lw=1.8, transform=ax.transAxes)
    note(ax, "约去共同因子", .34, color=RED)
    math(ax, r"=\frac{n!}{r!\,t!\,s!}", .20)
    save(fig, "公式03_组合数约分.png")


def formula03_coefficient():
    fig, ax = new_canvas(3.2)
    math(ax, r"\frac{n!}{r!t!s!}=\frac{3!}{2!1!0!}", .67)
    math(ax, r"=\frac{3\times2\times1}{(2\times1)\times1\times1}=3", .33)
    save(fig, "公式03_系数.png")


def formula05_mean_square():
    fig, ax = new_canvas(4.2)
    math(ax, r"E(Z_i)=1\cdot q_+ + 0\cdot q_0 + (-1)\cdot q_-", .73)
    math(ax, r"=q_+-q_-=\theta", .57)
    math(ax, r"E(Z_i^2)=1^2\cdot q_+ + 0^2\cdot q_0 + (-1)^2\cdot q_-", .34)
    math(ax, r"=q_++q_-", .16)
    save(fig, "公式05.png")


def formula05_substitution():
    fig, ax = new_canvas(3.5)
    math(ax, r"\mathrm{Var}(Z_i)=E(Z_i^2)-[E(Z_i)]^2", .74)
    math(ax, r"=(q_++q_-)-(q_+-q_-)^2", .49)
    math(ax, r"=q_++q_--\theta^2", .24)
    save(fig, "公式05_方差代入.png")


def formula06_average():
    fig, ax = new_canvas(9.5)
    math(ax, r"D=\frac{1}{n}\sum_{i=1}^{n}Z_i", .855)
    math(ax, r"E(D)=\frac{1}{n}\sum_{i=1}^{n}E(Z_i)=\theta", .715)
    math(ax, r"\mathrm{Var}(D)=\frac{1}{n^2}\mathrm{Var}\!\left(\sum_{i=1}^{n}Z_i\right)", .56)
    math(ax, r"=\frac{1}{n^2}\sum_{i=1}^{n}\mathrm{Var}(Z_i)", .405)
    note(ax, "不同人独立", .315)
    math(ax, r"=\frac{1}{n^2}\cdot n(q_++q_--\theta^2)", .245)
    math(ax, r"=\frac{q_++q_--\theta^2}{n}", .135)
    math(ax, r"\mathrm{SE}(D)=\sqrt{\frac{q_++q_--\theta^2}{n}}", .02)
    save(fig, "公式06.png")


if __name__ == "__main__":
    formula03_cancellation(); formula03_coefficient(); formula05_mean_square()
    formula05_substitution(); formula06_average()
