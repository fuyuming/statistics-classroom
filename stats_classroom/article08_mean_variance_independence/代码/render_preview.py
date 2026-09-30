#!/usr/bin/env python3
"""公众号 Markdown 子集 -> 手机预览。仅Python标准库。
显式URL自动可点击；代码保持原样；未知LaTeX报错，不猜译。
"""
import argparse
import html
import json
import os
import re
from pathlib import Path
from urllib.parse import quote, unquote, urlparse

CSS = '*{box-sizing:border-box}body{margin:0;background:#f1f4f3;color:#27343a;font-family:-apple-system,BlinkMacSystemFont,"PingFang SC","Microsoft YaHei",sans-serif}main{max-width:720px;margin:24px auto;padding:32px 36px 56px;background:#fff;border-top:5px solid #167d80}.brand{display:flex;align-items:center;gap:12px;color:#167d80;font-size:14px;letter-spacing:.04em}.brand img{width:48px;height:48px;border-radius:50%}h1{font-size:30px;line-height:1.45;color:#142e41;margin:26px 0 20px}h2{font-size:22px;line-height:1.5;color:#163d48;border-left:4px solid #198c8b;padding-left:13px;margin:42px 0 20px}h3{font-size:19px;line-height:1.55;color:#163d48;margin:30px 0 14px}p{font-size:17px;line-height:1.95;margin:19px 0;text-align:justify}strong{color:#15696a;font-weight:650}blockquote{margin:20px 0 30px;padding:10px 18px;background:#f0f7f6;border-left:2px solid #6da6a0}blockquote p{font-size:14px;line-height:1.9;margin:8px 0;color:#53666c}a{color:#167d80;text-decoration:none;border-bottom:1px solid #bbd8d4}code{font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:15px;background:#f2f7f6;color:#15696a;padding:1px 5px;border-radius:4px}pre{background:#f7faf9;border:1px solid #e4ecea;border-radius:6px;padding:12px 14px;overflow-x:auto}pre code{background:none;padding:0;font-size:14px;line-height:1.7}div.formula{background:#f7faf9;border-left:3px solid #6da6a0;padding:12px 16px;margin:20px 0;font-family:ui-monospace,SFMono-Regular,Menlo,monospace;font-size:17px;color:#15696a;text-align:center;line-height:1.9}ul,ol{font-size:17px;line-height:1.95;margin:18px 0;padding-left:26px}li{margin:8px 0}table{width:100%;border-collapse:collapse;margin:24px 0;font-size:15px}th,td{border-bottom:1px solid #e4ecea;padding:9px 10px;text-align:left;line-height:1.7}th{background:#f0f7f6;color:#15696a;font-weight:650}footer{border-top:1px solid #dfe8e6;margin-top:36px;padding-top:16px;font-size:12px;color:#758787;line-height:1.8}@media(max-width:600px){body{background:#fff}main{margin:0;padding:24px 22px 40px;border-top-width:4px}h1{font-size:27px}h2{font-size:21px}h3{font-size:18px}p,ul,ol{font-size:16px;line-height:1.95}table{font-size:14px}}main > p img{display:block;width:100%;height:auto;margin:24px 0 8px;border:1px solid #e4ecea}main > p:has(img){margin:25px -12px 4px}main > p:has(img) a{border:0;display:block}main > p > em:only-child{display:block;font-size:13px;line-height:1.8;font-style:normal;color:#657a7d;text-align:left;margin:0 0 28px}@media(max-width:600px){main > p:has(img){margin-left:-12px;margin-right:-12px}}a{overflow-wrap:anywhere}table{display:block;overflow-x:auto}blockquote p{font-size:16px}hr{border:0;border-top:1px solid #dfe8e6;margin:30px 0}'

# 本篇包含较长的代码目录名，窄屏允许换行，避免撑宽整页。
CSS += 'p,li,blockquote,h1,h2,h3{overflow-wrap:anywhere}main{min-width:0}'

# 只转换含义明确且经核对的少量公式；其他公式改为Unicode或公式图片。
FORMULAS = {
    r'n-1': 'n − 1',
    r'\mathrm{CV}=\frac{s}{\bar{x}}\times100\%': 'CV = s / x̄ × 100%',
    r'\mathrm{GSD}=\exp\left(s_{\ln x}\right)': 'GSD = exp(s_ln x)',
}

def formula_text(source):
    key = ''.join(source.split())
    if key in FORMULAS:
        return FORMULAS[key]
    # 已使用Unicode的公式不需要降级转换。
    if not re.search(r'[\\{}_^]', source):
        return source.strip()
    raise ValueError('未支持的LaTeX公式，请人工改成准确的Unicode或公式图片：' + source[:140])


def encode_url(url):
    parsed = urlparse(url)
    if parsed.scheme and parsed.scheme not in ('http', 'https', 'mailto'):
        raise ValueError('不支持的链接协议：' + parsed.scheme)
    return url if parsed.scheme else quote(unquote(url), safe='/#?=&')


def inline(source):
    """先保护代码、图片和链接，再处理正文URL，避免把代码或href二次转换。"""
    fragments = []

    def hold(value):
        index = len(fragments)
        fragments.append(value)
        return f'\x00{index}\x00'

    # 反引号中的**、$和网址必须保持代码原样。
    text = re.sub(r'`([^`]+)`', lambda m: hold('<code>' + html.escape(m[1]) + '</code>'), source)
    text = re.sub(
        r'!\[([^\]]*)\]\(([^)]+)\)',
        lambda m: hold('<img src="' + html.escape(encode_url(m[2]), quote=True) +
                       '" alt="' + html.escape(m[1], quote=True) + '">'), text)
    text = re.sub(
        r'\[([^\]]+)\]\(([^)]+)\)',
        lambda m: hold('<a href="' + html.escape(encode_url(m[2]), quote=True) +
                       '" target="_blank" rel="noopener">' + html.escape(m[1]) + '</a>'), text)
    text = re.sub(r'\$([^$]+)\$', lambda m: hold('<span class="formula-inline">' +
                   html.escape(formula_text(m[1])) + '</span>'), text)

    def explicit_link(match):
        url = match[0].rstrip('.,;:!?')
        trailing = match[0][len(url):]
        return hold('<a href="' + html.escape(url, quote=True) +
                    '" target="_blank" rel="noopener">' + html.escape(url) + '</a>') + trailing

    text = re.sub(r'https?://[A-Za-z0-9./_?=&%#:+~()@!-]+', explicit_link, text)
    # 第06篇使用sup表达幂次；只允许这几个无属性展示标签。
    text = re.sub(r'</?(?:sup|sub)>|<br\s*/?>', lambda m: hold(m[0]), text)
    text = html.escape(text)
    text = re.sub(r'\*\*([^*]+)\*\*', r'<strong>\1</strong>', text)
    text = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'<em>\1</em>', text)
    # 支持 [![图片](图.png)](图.png) 这类嵌套展示。
    for _ in range(len(fragments) + 1):
        if '\x00' not in text:
            break
        text = re.sub(r'\x00(\d+)\x00', lambda m: fragments[int(m[1])], text)
    return text


def convert(source, output, footer='', avatar=None, brand='明哥的微生物世界 · 统计课堂'):
    source = Path(source).resolve()
    output = Path(output).resolve()
    lines = source.read_text(encoding='utf-8').splitlines()
    body = []
    title = ''
    i = 0
    while i < len(lines):
        st = lines[i].strip()
        if not st:
            i += 1
            continue
        if st.startswith('```'):
            code = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith('```'):
                code.append(lines[i])
                i += 1
            if i == len(lines):
                raise ValueError('代码围栏未关闭')
            body.append('<pre><code>' + html.escape('\n'.join(code)) + '</code></pre>')
            i += 1
            continue
        if st.startswith('$$'):
            if st != '$$' and st.endswith('$$'):
                formula = st[2:-2]
                i += 1
            else:
                i += 1
                parts = []
                while i < len(lines) and lines[i].strip() != '$$':
                    parts.append(lines[i])
                    i += 1
                if i == len(lines):
                    raise ValueError('公式块未关闭')
                formula = ' '.join(parts)
                i += 1
            body.append('<div class="formula">' + html.escape(formula_text(formula)) + '</div>')
            continue
        if st.startswith('# ') and not title:
            title = st[2:]
            i += 1
            continue
        heading = re.match(r'^(#{2,3})\s+(.+)$', st)
        if heading:
            level = len(heading[1])
            body.append(f'<h{level}>' + inline(heading[2]) + f'</h{level}>')
            i += 1
            continue
        if st.startswith('>'):
            parts = []
            while i < len(lines) and lines[i].strip().startswith('>'):
                item = lines[i].strip()[1:].strip()
                if item:
                    parts.append('<p>' + inline(item) + '</p>')
                i += 1
            body.append('<blockquote>' + ''.join(parts) + '</blockquote>')
            continue
        if st == '---':
            body.append('<hr>')
            i += 1
            continue
        if st.startswith('|'):
            rows = []
            while i < len(lines) and lines[i].strip().startswith('|'):
                cells = re.split(r'(?<!\\)\|', lines[i].strip().strip('|'))
                cells = [c.strip().replace(r'\|', '|') for c in cells]
                if not set(''.join(cells)) <= set('-: '):
                    rows.append(cells)
                i += 1
            if rows:
                table = '<table><thead><tr>' + ''.join('<th>'+inline(c)+'</th>' for c in rows[0]) + '</tr></thead><tbody>'
                for row in rows[1:]:
                    if len(row) != len(rows[0]):
                        raise ValueError('表格列数不一致，竖线请用\\|转义')
                    table += '<tr>' + ''.join('<td>'+inline(c)+'</td>' for c in row) + '</tr>'
                body.append(table + '</tbody></table>')
            continue
        list_match = re.match(r'^([-*]|\d+\.)\s+(.+)', st)
        if list_match:
            ordered = list_match[1][0].isdigit()
            tag = 'ol' if ordered else 'ul'
            pattern = r'^\d+\.\s+(.+)' if ordered else r'^[-*]\s+(.+)'
            items = []
            while i < len(lines):
                item = re.match(pattern, lines[i].strip())
                if not item:
                    break
                items.append('<li>' + inline(item[1]) + '</li>')
                i += 1
            body.append('<'+tag+'>'+''.join(items)+'</'+tag+'>')
            continue
        body.append('<p>' + inline(st) + '</p>')
        i += 1
    if not title:
        raise ValueError('缺少一级标题')

    # 自动头像路径仅在真实存在时使用；可显式指定其他头像。
    avatar_path = Path(avatar).expanduser().resolve() if avatar else source.parent.parent / '账号视觉' / '明哥的微生物世界_头像_v1.png'
    avatar_html = ''
    if avatar_path.exists():
        relative = os.path.relpath(avatar_path, output.parent)
        avatar_html = '<img src="' + quote(relative) + '" alt="账号头像">'
    elif avatar:
        raise FileNotFoundError(avatar_path)
    document = ('<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">'
                '<meta name="viewport" content="width=device-width,initial-scale=1">'
                '<title>' + html.escape(title) + '</title><style>' + CSS + '</style></head><body><main>'
                '<div class="brand">' + avatar_html + html.escape(brand) + '</div><h1>' + inline(title) + '</h1>' +
                '\n'.join(body) + '<footer>' + html.escape(footer) + '</footer></main></body></html>')

    # 图片相对路径以稿件为起点；输出移到其他目录时，相应重写引用。
    image_count = 0
    def image_path(match):
        nonlocal image_count
        url = html.unescape(match[1])
        if urlparse(url).scheme in ('http', 'https'):
            image_count += 1
            return match[0]
        path = (source.parent / unquote(url)).resolve()
        # 头像已经按输出目录生成相对路径。
        if 'alt="账号头像"' in match[0]:
            path = (output.parent / unquote(url)).resolve()
        else:
            image_count += 1
        if not path.is_file():
            raise FileNotFoundError('图片引用不存在：' + str(path))
        rewritten = quote(os.path.relpath(path, output.parent))
        return match[0].replace(match[1], rewritten, 1)
    document = re.sub(r'<img src="([^"]+)"[^>]*>', image_path, document)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(document, encoding='utf-8')
    return {'title': title, 'content_images': image_count, 'source': str(source),
            'output': str(output), 'characters_in_markdown': len(source.read_text(encoding='utf-8'))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('output')
    parser.add_argument('--footer', default='本地预览 · 尚未发布')
    parser.add_argument('--avatar')
    parser.add_argument('--brand', default='明哥的微生物世界 · 统计课堂')
    args = parser.parse_args()
    result = convert(args.source, args.output, args.footer, args.avatar, args.brand)
    print(json.dumps(result, ensure_ascii=False, indent=2))
