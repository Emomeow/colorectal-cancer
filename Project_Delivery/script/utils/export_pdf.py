#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
通用报告 PDF 生成脚本
用法: python export_pdf.py <input.md> <output.pdf>
"""
import sys, io, re, base64, subprocess, tempfile, argparse
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')
import markdown

BASE_DIR = Path('d:/changdao')
EDGE     = Path('C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe')

parser = argparse.ArgumentParser()
parser.add_argument('md_file')
parser.add_argument('out_pdf')
args = parser.parse_args()

MD_FILE = Path(args.md_file)
OUT_PDF = Path(args.out_pdf)

md_text = MD_FILE.read_text(encoding='utf-8')

def embed_img(m):
    alt, path = m.group(1), m.group(2)
    if path.startswith('http') or path.startswith('data:'):
        return m.group(0)
    img_path = (BASE_DIR / path).resolve()
    if not img_path.exists():
        print(f'  [WARN] missing: {img_path}')
        return f'<p><em>[missing: {alt}]</em></p>'
    ext  = img_path.suffix.lower().lstrip('.')
    mime = {'jpg': 'jpeg', 'jpeg': 'jpeg', 'png': 'png'}.get(ext, 'png')
    b64  = base64.b64encode(img_path.read_bytes()).decode('ascii')
    return f'<img alt="{alt}" src="data:image/{mime};base64,{b64}" style="max-width:100%">'

md_text = re.sub(r'!\[([^\]]*)\]\(([^)]+)\)', embed_img, md_text)
html_body = markdown.markdown(md_text, extensions=['tables', 'fenced_code'])

html = f'''<!DOCTYPE html>
<html><head><meta charset="utf-8">
<style>
@page {{
  size: A4;
  margin: 15mm 12mm 15mm 12mm;
}}
body {{
  font-family: Arial, "Microsoft YaHei", sans-serif;
  font-size: 11pt;
  line-height: 1.65;
  color: #222;
  margin: 0;
  padding: 0 8px;
}}
h1 {{ font-size: 18pt; color: #1a1a2e; border-bottom: 2px solid #1a1a2e; padding-bottom: 5pt; margin-top: 0; }}
h2 {{ font-size: 14pt; color: #16213e; border-bottom: 1px solid #ccc; padding-bottom: 3pt; margin-top: 24pt; }}
h3 {{ font-size: 12pt; color: #0f3460; margin-top: 14pt; }}
h4 {{ font-size: 11pt; color: #555; margin-top: 8pt; }}
table {{ border-collapse: collapse; width: 100%; margin: 8pt 0; font-size: 9.5pt; }}
th {{ background: #1a1a2e; color: white; padding: 5pt 8pt; text-align: left; }}
td {{ border: 1px solid #ddd; padding: 4pt 8pt; }}
tr:nth-child(even) td {{ background: #f7f7f7; }}
img {{ max-width: 100%; height: auto; display: block; margin: 10pt auto; border: 1px solid #e8e8e8; }}
blockquote {{ border-left: 4px solid #4DBBD5; margin: 6pt 0; padding: 4pt 12pt; background: #f0f8ff; color: #444; }}
code {{ background: #f4f4f4; padding: 1pt 4pt; border-radius: 3px; font-size: 9pt; font-family: Consolas, monospace; }}
pre code {{ display: block; padding: 8pt; overflow-x: auto; }}
hr {{ border: none; border-top: 1px solid #ddd; margin: 14pt 0; }}
p {{ margin: 4pt 0; }}
ul, ol {{ margin: 4pt 0; padding-left: 18pt; }}
li {{ margin: 2pt 0; }}
strong {{ color: #1a1a2e; }}
@media print {{
  img {{ page-break-inside: avoid; max-height: 200mm; }}
  table {{ page-break-inside: avoid; }}
}}
</style>
</head><body>
{html_body}
</body></html>'''

with tempfile.NamedTemporaryFile(suffix='.html', delete=False, mode='w', encoding='utf-8') as f:
    f.write(html)
    tmp = f.name

print(f'Temp HTML: {tmp}')
print('Converting via Edge headless ...')

cmd = [
    str(EDGE),
    '--headless=new',
    '--disable-gpu',
    '--no-sandbox',
    '--no-pdf-header-footer',
    '--print-to-pdf-no-header',
    f'--print-to-pdf={OUT_PDF.resolve()}',
    tmp,
]
subprocess.run(cmd, capture_output=True, timeout=420)
Path(tmp).unlink(missing_ok=True)

if OUT_PDF.exists() and OUT_PDF.stat().st_size > 1000:
    print(f'Saved: {OUT_PDF}  ({OUT_PDF.stat().st_size/1e6:.1f} MB)')
else:
    print('ERROR: PDF not created')
