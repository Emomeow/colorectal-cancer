#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate report_phase1.pdf from report_phase1.md using reportlab.
Supports Chinese text via Microsoft YaHei font.
"""

import re
import os
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.lib import colors
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    HRFlowable, KeepTogether
)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# ── Register fonts ─────────────────────────────────────────────────────────────
pdfmetrics.registerFont(TTFont('MSYaHei', 'C:/Windows/Fonts/msyh.ttc', subfontIndex=0))
pdfmetrics.registerFont(TTFont('MSYaHeiBold', 'C:/Windows/Fonts/msyhbd.ttc', subfontIndex=0))
pdfmetrics.registerFont(TTFont('Arial', 'C:/Windows/Fonts/Arial.ttf'))
try:
    pdfmetrics.registerFont(TTFont('ArialBold', 'C:/Windows/Fonts/Arialbd.ttf'))
except Exception:
    pdfmetrics.registerFont(TTFont('ArialBold', 'C:/Windows/Fonts/Arial.ttf'))

MAIN_FONT = 'MSYaHei'
BOLD_FONT = 'MSYaHeiBold'

# ── Styles ─────────────────────────────────────────────────────────────────────
styles = getSampleStyleSheet()

def make_style(name, parent='Normal', font=None, **kwargs):
    fn = font or MAIN_FONT
    return ParagraphStyle(name, parent=styles[parent], fontName=fn, **kwargs)

style_h1    = make_style('H1',    font=BOLD_FONT,  fontSize=16, spaceAfter=6, leading=22,
                         textColor=colors.HexColor('#1a1a2e'))
style_h2    = make_style('H2',    font=BOLD_FONT,  fontSize=13, spaceAfter=4, spaceBefore=14,
                         leading=18, textColor=colors.HexColor('#16213e'))
style_h3    = make_style('H3',    font=BOLD_FONT,  fontSize=11, spaceAfter=3, spaceBefore=8,
                         leading=16, textColor=colors.HexColor('#0f3460'))
style_body  = make_style('Body',  fontSize=9.5, spaceAfter=3, leading=14)
style_quote = make_style('Quote', fontSize=9,   spaceAfter=4, leading=13,
                         leftIndent=14, textColor=colors.HexColor('#555555'))
style_bullet = make_style('Bullet', fontSize=9.5, spaceAfter=2, leading=13,
                           leftIndent=12, firstLineIndent=-10)
style_code  = make_style('Code',  font='Arial', fontSize=8.5, spaceAfter=3, leading=12,
                         leftIndent=12)

TABLE_STYLE = TableStyle([
    ('BACKGROUND',  (0, 0), (-1, 0),  colors.HexColor('#2c3e7a')),
    ('TEXTCOLOR',   (0, 0), (-1, 0),  colors.white),
    ('FONTNAME',    (0, 0), (-1, 0),  BOLD_FONT),
    ('FONTSIZE',    (0, 0), (-1, 0),  8.5),
    ('FONTNAME',    (0, 1), (-1, -1), MAIN_FONT),
    ('FONTSIZE',    (0, 1), (-1, -1), 8.5),
    ('ROWBACKGROUNDS', (0, 1), (-1, -1),
     [colors.white, colors.HexColor('#f0f4ff')]),
    ('GRID',        (0, 0), (-1, -1), 0.4, colors.HexColor('#cccccc')),
    ('VALIGN',      (0, 0), (-1, -1), 'MIDDLE'),
    ('TOPPADDING',  (0, 0), (-1, -1), 4),
    ('BOTTOMPADDING', (0, 0), (-1, -1), 4),
    ('LEFTPADDING', (0, 0), (-1, -1), 6),
    ('RIGHTPADDING', (0, 0), (-1, -1), 6),
    ('WORDWRAP',    (0, 0), (-1, -1), True),
])

# ── Markdown parser → ReportLab flowables ─────────────────────────────────────
def inline_fmt(text):
    """Convert inline markdown (bold, code) to ReportLab XML."""
    # Bold: **text** or __text__
    text = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', text)
    text = re.sub(r'__(.+?)__',     r'<b>\1</b>', text)
    # Italic: *text* or _text_
    text = re.sub(r'\*(.+?)\*', r'<i>\1</i>', text)
    # Code: `text`
    text = re.sub(r'`(.+?)`',
                  r'<font name="Arial" size="8.5" color="#c7254e">\1</font>', text)
    # Escape bare & < > that aren't part of tags
    # (simple approach: already replaced bold/italic, leftover < > need escaping)
    # Only escape if not inside an existing tag
    def escape_bare(m):
        return m.group(0).replace('&', '&amp;') if '&amp;' not in m.group(0) else m.group(0)
    return text

def md_to_flowables(md_text):
    flowables = []
    lines = md_text.split('\n')
    table_rows = []
    in_table = False
    i = 0

    while i < len(lines):
        line = lines[i]

        # Table row detection
        if line.strip().startswith('|'):
            if not in_table:
                in_table = True
                table_rows = []
            # Skip separator rows (---|---|---)
            if re.match(r'^\|[\s\-|:]+\|$', line.strip()):
                i += 1
                continue
            cells = [c.strip() for c in line.strip().strip('|').split('|')]
            table_rows.append(cells)
            i += 1
            continue
        else:
            if in_table and table_rows:
                # Flush table
                flowables.append(_make_table(table_rows))
                flowables.append(Spacer(1, 6))
                table_rows = []
                in_table = False

        # Headings
        if line.startswith('### '):
            flowables.append(Spacer(1, 4))
            flowables.append(Paragraph(inline_fmt(line[4:]), style_h3))
        elif line.startswith('## '):
            flowables.append(Spacer(1, 8))
            flowables.append(HRFlowable(width='100%', thickness=1.5,
                                         color=colors.HexColor('#2c3e7a'), spaceAfter=4))
            flowables.append(Paragraph(inline_fmt(line[3:]), style_h2))
        elif line.startswith('# '):
            flowables.append(Paragraph(inline_fmt(line[2:]), style_h1))
            flowables.append(HRFlowable(width='100%', thickness=2,
                                         color=colors.HexColor('#1a1a2e'), spaceAfter=6))
        # Blockquote
        elif line.startswith('> '):
            text = inline_fmt(line[2:])
            flowables.append(Paragraph(text, style_quote))
        # Horizontal rule
        elif re.match(r'^---+\s*$', line):
            flowables.append(Spacer(1, 4))
            flowables.append(HRFlowable(width='100%', thickness=0.5,
                                         color=colors.HexColor('#cccccc'), spaceAfter=4))
        # Bullet list
        elif line.startswith('- ') or line.startswith('* '):
            text = '• ' + inline_fmt(line[2:])
            flowables.append(Paragraph(text, style_bullet))
        # Numbered list
        elif re.match(r'^\d+\. ', line):
            text = inline_fmt(line)
            flowables.append(Paragraph(text, style_bullet))
        # Image (skip — figures not embedded)
        elif line.startswith('!['):
            alt = re.search(r'!\[(.+?)\]', line)
            alt_text = alt.group(1) if alt else 'Figure'
            flowables.append(Paragraph(f'<i>[图片：{alt_text}]</i>', style_body))
        # Empty line
        elif line.strip() == '':
            if flowables and not isinstance(flowables[-1], Spacer):
                flowables.append(Spacer(1, 4))
        # Bold-only line (treat as mini-heading)
        elif re.match(r'^\*\*.+\*\*$', line.strip()):
            flowables.append(Paragraph(inline_fmt(line.strip()), style_h3))
        # Normal paragraph
        else:
            if line.strip():
                flowables.append(Paragraph(inline_fmt(line.strip()), style_body))

        i += 1

    # Flush remaining table
    if in_table and table_rows:
        flowables.append(_make_table(table_rows))

    return flowables


def _make_table(rows):
    """Build a ReportLab Table from list-of-lists."""
    if not rows:
        return Spacer(1, 0)

    n_cols = max(len(r) for r in rows)
    # Pad rows to equal width
    padded = [r + [''] * (n_cols - len(r)) for r in rows]

    # Wrap cells in Paragraph for word wrap + Chinese support
    wrapped = []
    for ri, row in enumerate(padded):
        style = make_style(f'tcell_{ri}',
                           font=BOLD_FONT if ri == 0 else MAIN_FONT,
                           fontSize=8.5, leading=12,
                           textColor=colors.white if ri == 0 else colors.black)
        wrapped.append([Paragraph(inline_fmt(str(c)), style) for c in row])

    page_width = A4[0] - 4 * cm
    col_width = page_width / n_cols
    tbl = Table(wrapped, colWidths=[col_width] * n_cols, repeatRows=1)
    tbl.setStyle(TABLE_STYLE)
    return tbl


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    md_path = 'report_phase1.md'
    pdf_path = 'report_phase1.pdf'

    print(f"Reading {md_path}...")
    with open(md_path, encoding='utf-8') as f:
        md_text = f.read()

    print("Converting to PDF flowables...")
    flowables = md_to_flowables(md_text)

    print(f"Writing {pdf_path}...")
    doc = SimpleDocTemplate(
        pdf_path,
        pagesize=A4,
        leftMargin=2*cm, rightMargin=2*cm,
        topMargin=2*cm, bottomMargin=2*cm,
        title='阶段一报告 — CRC/UC 免疫微环境 scRNA-seq 整合分析',
        author='changdao project',
    )
    doc.build(flowables)
    size_mb = os.path.getsize(pdf_path) / 1024**2
    print(f"Done: {pdf_path}  ({size_mb:.1f} MB)")


if __name__ == '__main__':
    main()
