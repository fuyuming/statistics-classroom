from pathlib import Path
from md2html_preview import convert
root=Path(__file__).resolve().parent.parent
out=root/'一条直线_手机预览_v5.html'
convert(root/'一条直线_公众号稿_v5.md',out,footer='Statistical Rethinking 精读 03 · v5 同一模型两种计算')
(root/'一条直线_手机预览.html').write_text(out.read_text())
