from pathlib import Path
from md2html_preview import convert
root=Path(__file__).resolve().parent.parent
out=root/'一条直线_手机预览_v4.html'
convert(root/'一条直线_公众号稿_v4.md',out,footer='Statistical Rethinking 精读 03 · v4 配图版')
(root/'一条直线_手机预览.html').write_text(out.read_text())
