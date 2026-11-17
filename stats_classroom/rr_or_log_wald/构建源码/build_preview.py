"""从唯一源稿生成普通与单文件手机预览；仅Python标准库。
运行：python 构建源码/build_preview.py。重新生成配图时先运行build_assets.py。
"""
from pathlib import Path
import base64
import re
from urllib.parse import unquote
import md2html_preview

root = Path(__file__).resolve().parents[1]
source = root / "RR与OR_公众号稿.md"
output = root / "手机预览.html"
md2html_preview.convert(source, output,
    brand="明哥的微生物世界 · 统计课堂（篇号待定）")
content = output.read_text(encoding="utf-8")
# 使长标题各行长度均衡，避免一个字单独占行。
content = content.replace("</style>", "h1{text-wrap:balance}p{text-align:left}@media(max-width:600px){h1{font-size:25px}}</style>")
output.write_text(content, encoding="utf-8")

def embed(match):
    name = unquote(match.group(1))
    path = root / name
    blob = base64.b64encode(path.read_bytes()).decode("ascii")
    return f'src="data:image/png;base64,{blob}"'

single = re.sub(r'src="([^\"]+\.png)"', embed, content)
(root / "手机预览_单文件.html").write_text(single, encoding="utf-8")
print("普通预览与单文件预览已从同一源稿生成。")
