"""从当前Markdown重建离线手机预览。只需latex2mathml，已有HTML可直接打开。"""
from pathlib import Path
import subprocess, json
root=Path(__file__).resolve().parent.parent
import importlib.util, re, html, sys, struct
from latex2mathml.converter import convert as latex_to_mathml
import xml.etree.ElementTree as ET
spec=importlib.util.spec_from_file_location('preview_base',Path(__file__).with_name('preview_base.py'))
base=importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
original_inline=base.inline
math_inventory=[]
def typeset_inline(source):
    saved=[]
    def protect(match):
        if match[0].startswith('`'):
            return match[0]
        tex=match[0][1:-1]
        key=f'MATHSAFEPLACEHOLDER{len(saved)}END'
        math_root=ET.fromstring(latex_to_mathml(tex,display='inline'))
        for node in math_root.iter():
            if node.tag.endswith('mover') and len(node)==2 and node[1].text=='^':
                node.set('accent','true')
                node[1].set('stretchy','true')
                node[1].text='\u02c6'
        mml=ET.tostring(math_root,encoding='unicode').replace('ns0:','').replace(':ns0','')
        saved.append('<span class="inline-math" data-tex="'+html.escape(tex,quote=True)+'">'+mml+'</span>')
        math_inventory.append(tex)
        return key
    protected=re.sub(r'`[^`]+`|\$[^$\n]+\$',protect,source)
    rendered=original_inline(protected)
    for i,value in enumerate(saved):
        rendered=rendered.replace(f'MATHSAFEPLACEHOLDER{i}END',value)
    rendered=re.sub(r'</math></span>(分(?:<sup>[24]</sup>)?)',r'</math>\1</span>',rendered)
    rendered=re.sub(r'</span>([，。、；：？！])',r'\1</span>',rendered)
    return rendered
base.inline=typeset_inline
base.convert(root/'双向有序条件模拟_公众号稿.md',root/'手机预览.html',brand='明哥的微生物世界 · 统计课堂 41')
# 针对本篇六列表和公式统一字号，沿用白底、深色正文与亮蓝标题。
p=root/'手机预览.html';document=p.read_text();document=document.replace('<p>明哥的微生物世界 · 统计课堂 41</p>', '', 1);document=document.replace('</style>','.inline-math{display:inline-block;white-space:nowrap;font-weight:normal;max-width:100%;vertical-align:baseline}.inline-math math{font-family:"STIX Two Math","Cambria Math",serif;font-size:1.04em}strong .inline-math{font-weight:bold} @media(max-width:600px){table{display:table;table-layout:fixed;font-size:12px;width:100%}th,td{padding:7px 3px;text-align:center;overflow-wrap:anywhere}th:first-child,td:first-child{width:22%}pre{white-space:pre-wrap;overflow-wrap:anywhere}pre code{font-size:12px}p{text-align:left}code{overflow-wrap:anywhere;word-break:break-word;white-space:normal}pre code{white-space:pre-wrap;word-break:normal}}\n</style>');p.write_text(document)

# 数学图按同一显示字号定宽；PNG的高分辨率像素不直接作为CSS像素。
layout = json.loads(Path(__file__).with_name('公式样式.json').read_text())
from urllib.parse import unquote
def size_formula(match):
    tag = match[0]
    src = re.search(r'\bsrc="([^"]+)"', tag)
    if not src:
        return tag
    image_path = root / unquote(html.unescape(src[1]))
    if not image_path.name.startswith('公式') or image_path.suffix.lower() != '.png':
        return tag
    header = image_path.read_bytes()[:24]
    if header[:8] != b'\x89PNG\r\n\x1a\n':
        raise ValueError(f'公式图不是有效PNG：{image_path.name}')
    pixels_w, _ = struct.unpack('>II', header[16:24])
    display_w = pixels_w / layout['device_scale_factor'] * layout['display_scale']
    display_font = layout['font_size_css'] * layout['display_scale']
    attrs = (f' class="formula-image" data-math-font-px="{display_font:g}"'
             f' style="display:block;width:{display_w:g}px;max-width:100%;height:auto;'
             'margin-left:auto;margin-right:auto;border:0"')
    return tag[:-1] + attrs + '>'
document = re.sub(r'<img\b[^>]*>', size_formula, document)
p.write_text(document)

print(root/'手机预览.html')
