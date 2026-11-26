"""从公众号稿.md唯一正文源生成离线HTML与公式图；只在本地运行。
依赖：matplotlib, Pillow, lxml；使用同目录提供的预览渲染器。
"""
from pathlib import Path
import re
import os
import json
import base64
import importlib.util
from urllib.parse import unquote
from PIL import Image
from lxml import html, etree
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'公众号稿.md'
FORMULAS=ROOT/'文章配图/公式'
FORMULAS.mkdir(exist_ok=True)
WORK=ROOT.parents[1]/'work/residual_preview_build'
WORK.mkdir(parents=True,exist_ok=True)
ASSET_REF=Path(os.path.relpath(ROOT/"文章配图",WORK)).as_posix()
manifest=[]

def formula(match):
    latex=' '.join(match.group(1).split())
    dest=FORMULAS/f'公式_{len(manifest)+1:02d}.png'
    math_to_image('$'+latex+'$',str(dest),prop=FontProperties(size=15),dpi=288,color='#263544')
    w,h=Image.open(dest).size
    # 全部公式用相同15pt；手机CSS按96dpi还原，不逐式缩字号。
    width=w/3
    if width>340:
        raise ValueError(f'公式宽度{width:.1f}超出手机正文，请拆行：{latex}')
    manifest.append({'file':dest.name,'latex':latex,'css_width':width,'css_height':h/3})
    return f'![公式：{latex}]({ASSET_REF}/公式/{dest.name})'

text=SOURCE.read_text(encoding='utf-8')
# 页首品牌由渲染器统一输出；源稿保留系列行，派生稿避免重复显示。
text=re.sub(r'^明哥的微生物世界 · 统计课堂\s+\d+\s*\n+', '', text, count=1)
text=re.sub(r'(!\[[^\]]*\]\()文章配图/',lambda m: m.group(1)+ASSET_REF+'/',text)
text=re.sub(r'\$\$\s*\n(.*?)\n\$\$',formula,text,flags=re.S)
for old_formula in FORMULAS.glob('公式_*.png'):
    if old_formula.name not in {item['file'] for item in manifest}:
        old_formula.unlink()
derived=WORK/'渲染源.md'
derived.write_text(text,encoding='utf-8')
renderer_path=Path(__file__).resolve().parent/'预览渲染器.py'
spec=importlib.util.spec_from_file_location('renderer',renderer_path)
renderer=importlib.util.module_from_spec(spec)
spec.loader.exec_module(renderer)
dest=ROOT/'手机预览.html'
renderer.convert(derived,dest,brand='明哥的微生物世界 · 统计课堂 43')
doc=html.fromstring(dest.read_text())
css=etree.SubElement(doc.find('head'),'style')
css.text='''
*{box-sizing:border-box}body{background:#eef2f6;color:#263544}
main{max-width:500px;margin:24px auto;padding:26px 20px 42px;background:white}
h1,h2,h3{color:#0000F5}h1{font-size:25px;line-height:1.5}h2{font-size:21px;line-height:1.6}h3{font-size:18px;line-height:1.6}
p,li{font-size:16px;line-height:1.95;overflow-wrap:anywhere;text-align:left}strong{color:#263544}
pre{white-space:pre-wrap;overflow-wrap:anywhere;padding:12px;font-size:13px}pre code{font-size:13px;line-height:1.65;font-variant-ligatures:none}
code{font-size:14px;font-variant-ligatures:none}table{display:table;width:100%;table-layout:auto;font-size:13px;line-height:1.7}
th,td{padding:8px 5px;vertical-align:middle;text-align:left;overflow-wrap:anywhere}th{color:#263544;background:#edf4ff}
em{font-size:13px;line-height:1.65;display:block;color:#607286;font-style:normal}
img{max-width:100%;height:auto}main p img{margin:4px auto;display:block;border:0}
sup,sub{font-size:75%;line-height:0;white-space:nowrap;overflow-wrap:normal}u{text-underline-offset:3px;text-decoration-thickness:1px}a{overflow-wrap:anywhere;color:#2358ac}.brand{font-size:13px;letter-spacing:0}
@media(max-width:600px){body{background:white}main{margin:0 auto;width:100%;padding:22px 20px 36px}}
'''
byname={r['file']:r for r in manifest}
for el in doc.xpath('//img'):
    path=(ROOT/unquote(el.get('src'))).resolve()
    if path.name in byname:
        item=byname[path.name]
        el.set('class','math-image')
        el.set('style',f'width:{item["css_width"]:.3f}px;max-width:100%;height:auto;')
    el.set('src','data:image/png;base64,'+base64.b64encode(path.read_bytes()).decode())
dest.write_text('<!doctype html>\n'+etree.tostring(doc,encoding='unicode',method='html'),encoding='utf-8')
(WORK/'公式清单.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'formulas':len(manifest),'largest_formula_css_width':max(x['css_width'] for x in manifest),'html':str(dest),'bytes':dest.stat().st_size},ensure_ascii=False))
