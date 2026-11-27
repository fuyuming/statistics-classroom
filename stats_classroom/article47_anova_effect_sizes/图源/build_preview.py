"""从唯一Markdown源稿重建本地HTML；仅用Python标准库。"""
from pathlib import Path
from render_md import convert
BASE=Path(__file__).resolve().parents[1]
out=BASE/'手机预览.html'
print(convert(BASE/'公众号稿.md',out,brand='明哥的微生物世界 · 统计课堂 47'))
s=out.read_text(encoding='utf-8')
# 源稿保留完整系列行；页面页首已有品牌栏，避免显示两次。
s=s.replace('<p>明哥的微生物世界 · 统计课堂 47</p>','',1)
style='''body{background:#fff;color:#242830}p{text-align:left}code{overflow-wrap:anywhere;word-break:break-word}main{border-top:0}h1{color:#242830}h2,h3{color:#0000F5}sup{font-size:.7em;line-height:0}table{display:table;table-layout:fixed;width:100%}td,th{overflow-wrap:anywhere}table:first-of-type th:first-child{width:18%}pre{max-width:100%;white-space:pre-wrap;overflow-wrap:anywhere}pre code{font-size:12px;line-height:1.7}main > p img{border:0}main > p:has(img){margin-left:0;margin-right:0}@media(max-width:600px){main{padding:22px 18px 40px}h1{font-size:25px}h2{font-size:21px}main > p:has(img){margin-left:0;margin-right:0}td,th{padding:8px 6px}table{font-size:14px}}'''
s=s.replace('</style>',style+'</style>')
out.write_text(s,encoding='utf-8')
