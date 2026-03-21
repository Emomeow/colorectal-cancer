#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convert report_phase1.md → HTML → PDF via Microsoft Edge headless"""
import sys, io, re, base64, subprocess, tempfile
from pathlib import Path

if sys.platform == 'win32':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

import markdown

BASE_DIR = Path('d:/changdao')
MD_FILE  = BASE_DIR / 'report_phase1.md'
OUT_PDF  = BASE_DIR / 'report_phase1.pdf'
EDGE     = Path('C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe')

# ── Embed images as base64 ─────────────────────────────────────
md_text = MD_FILE.read_text(encoding='utf-8')

def embed_img(m):
    alt, path = m.group(1), m.group(2)
    if path.startswith('http') or path.startswith('data:'):
        return m.group(0)
    img_path = (BASE_DIR / path).resolve()
    if not img_path.exists():
        print(f'  [WARN] missing: {img_path}')
        return f'<p><em>[图片缺失: {alt}]</em></p>'
    ext  = img_path.suffix.lower().lstrip('.')
    mime = {'jpg': 'jpeg', 'jpeg': 'jpeg', 'png': 'png'}.get(ext, 'png')
    b64  = base64.b64encode(img_path.read_bytes()).decode('ascii')
    return f'![{alt}](data:image/{mime};base64,{b64})'

md_text = re.sub(r'!\[([^\]]*)\]\(([^)]+)\)', embed_img, md_text)

# ── Markdown → HTML body ───────────────────────────────────────
body_html = markdown.markdown(
    md_text,
    extensions=['tables', 'toc', 'fenced_code']
)

# ── Full HTML with CSS (Edge renders system fonts natively) ────
full_html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8"/>
<style>
  * {{ box-sizing: border-box; }}
  body {{
    font-family: "Microsoft YaHei", "微软雅黑", "SimHei", "黑体", Arial, sans-serif;
    font-size: 11pt;
    line-height: 1.65;
    color: #222;
    margin: 0;
    padding: 0 8px;
  }}
  h1 {{
    font-size: 20pt;
    color: #1a1a2e;
    border-bottom: 2.5px solid #1a1a2e;
    padding-bottom: 6pt;
    margin-top: 0;
  }}
  h2 {{
    font-size: 14pt;
    color: #16213e;
    border-bottom: 1px solid #ccc;
    padding-bottom: 4pt;
    margin-top: 28pt;
    page-break-before: auto;
  }}
  h3 {{
    font-size: 12pt;
    color: #0f3460;
    margin-top: 16pt;
  }}
  h4 {{
    font-size: 11pt;
    color: #555;
    margin-top: 10pt;
  }}
  table {{
    border-collapse: collapse;
    width: 100%;
    margin: 10pt 0;
    font-size: 9.5pt;
  }}
  th {{
    background: #1a1a2e;
    color: white;
    padding: 5pt 8pt;
    text-align: left;
    font-weight: bold;
  }}
  td {{
    border: 1px solid #ddd;
    padding: 4pt 8pt;
  }}
  tr:nth-child(even) td {{ background: #f7f7f7; }}
  img {{
    max-width: 100%;
    height: auto;
    display: block;
    margin: 12pt auto;
    border: 1px solid #e8e8e8;
  }}
  blockquote {{
    border-left: 4px solid #4DBBD5;
    margin: 8pt 0;
    padding: 4pt 14pt;
    background: #f0f8ff;
    color: #444;
  }}
  code {{
    background: #f4f4f4;
    padding: 1pt 4pt;
    border-radius: 3px;
    font-size: 9pt;
    font-family: "Consolas", "Courier New", monospace;
  }}
  pre code {{
    display: block;
    padding: 8pt;
    overflow-x: auto;
  }}
  hr {{
    border: none;
    border-top: 1px solid #ddd;
    margin: 16pt 0;
  }}
  p {{ margin: 5pt 0; }}
  ul, ol {{ margin: 5pt 0; padding-left: 20pt; }}
  li {{ margin: 2pt 0; }}
  strong {{ color: #1a1a2e; }}
  @media print {{
    body {{ padding: 0; }}
    h2 {{ page-break-before: auto; }}
    img {{ page-break-inside: avoid; max-height: 220mm; }}
    table {{ page-break-inside: avoid; }}
  }}
</style>
</head>
<body>
{body_html}
</body>
</html>"""

# ── Write temp HTML and call Edge headless ─────────────────────
with tempfile.NamedTemporaryFile(suffix='.html', mode='w',
                                  encoding='utf-8', delete=False) as f:
    f.write(full_html)
    tmp_html = Path(f.name)

print(f'Temp HTML: {tmp_html}')
print('Converting via Edge headless ...')

cmd = [
    str(EDGE),
    '--headless=new',
    '--disable-gpu',
    '--no-sandbox',
    '--no-pdf-header-footer',
    f'--print-to-pdf={OUT_PDF}',
    f'--print-to-pdf-no-header',
    str(tmp_html),
]

result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
tmp_html.unlink(missing_ok=True)

if OUT_PDF.exists() and OUT_PDF.stat().st_size > 1000:
    print(f'Saved: {OUT_PDF}  ({OUT_PDF.stat().st_size / 1024 / 1024:.1f} MB)')
else:
    print('ERROR: PDF not created or too small')
    print(result.stdout[:500])
    print(result.stderr[:500])
