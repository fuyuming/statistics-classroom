"""从当前正文生成包含封面的预览；只依赖 Python 标准库。"""
from pathlib import Path
from urllib.parse import quote
from md2html_preview import convert

root = Path(__file__).resolve().parent.parent
output = root / '数路径_手机预览_v2.html'
convert(root / '数路径_公众号稿_v2.md', output,
        footer='Statistical Rethinking 精读 02 · 2026-10-02 修订版',
        avatar=root / '文章配图/账号头像.png',
        brand='明哥的微生物世界 · Statistical Rethinking 精读 02')
document = output.read_text(encoding='utf-8')
cover = quote('文章配图/00-公众号封面-精读02-横版_v2.png')
document = document.replace('</h1>', '</h1><p><img src="' + cover + '" alt="精读02新版封面"></p>', 1)
output.write_text(document, encoding='utf-8')
(root / '数路径_手机预览.html').write_text(document, encoding='utf-8')
print('已生成 v2 及当前入口预览。')
