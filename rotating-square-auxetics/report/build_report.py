"""Inline every fig/*.png of report_src.html as base64 -> one self-contained HTML file.

Run from this folder:  python build_report.py
Output: ../Rotating_Squares_Report.html   (KaTeX still loads from cdn.jsdelivr.net)
"""
import base64
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
src = (HERE / "report_src.html").read_text(encoding="utf-8")


def inline(m):
    data = base64.b64encode((HERE / m.group(1)).read_bytes()).decode("ascii")
    return f'src="data:image/png;base64,{data}"'


out = re.sub(r'src="(fig/[^"]+\.png)"', inline, src)
dst = HERE.parent / "Rotating_Squares_Report.html"
dst.write_text(out, encoding="utf-8")
print(f"wrote {dst}  ({dst.stat().st_size/1e6:.1f} MB, {src.count('fig/')} images inlined)")
