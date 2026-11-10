"""Convert report.html to GitHub-flavoured Markdown (report.md).

Images stay as relative links to fig/, maths uses GitHub's $...$ / $$...$$ syntax.
Run:  python to_markdown.py
"""
import html
import re
from pathlib import Path

from bs4 import BeautifulSoup, NavigableString, Tag

HERE = Path(__file__).resolve().parent


def inline(node):
    """Inline content -> Markdown text."""
    out = []
    for c in node.children:
        if isinstance(c, NavigableString):
            out.append(str(c))
        elif c.name in ("b", "strong"):
            t = inline(c).strip()
            out.append(f"**{t}**" if t else "")
        elif c.name in ("i", "em"):
            t = inline(c).strip()
            out.append(f"*{t}*" if t else "")
        elif c.name == "span" and "no" in c.get("class", []):
            out.append(inline(c).strip() + " ")
        elif c.name == "sub":
            out.append(f"<sub>{inline(c)}</sub>")
        elif c.name == "sup":
            out.append(f"<sup>{inline(c)}</sup>")
        elif c.name == "br":
            out.append(" ")
        elif c.name == "img":
            pass
        else:
            out.append(inline(c))
    s = "".join(out)
    s = re.sub(r"\\\((.+?)\\\)", lambda m: f"${m.group(1).strip()}$", s, flags=re.S)   # inline maths
    return re.sub(r"\s+", " ", s)


def table(t):
    rows = []
    for tr in t.find_all("tr"):
        cells = [inline(td).strip().replace("|", "\\|") for td in tr.find_all(["td", "th"])]
        rows.append(cells)
    if not rows:
        return ""
    n = max(len(r) for r in rows)
    rows = [r + [""] * (n - len(r)) for r in rows]
    md = ["| " + " | ".join(rows[0]) + " |", "|" + "---|" * n]
    md += ["| " + " | ".join(r) + " |" for r in rows[1:]]
    return "\n".join(md)


def figure(f):
    imgs = f.find_all("img")
    cap = f.find("figcaption")
    capt = inline(cap).strip() if cap else ""
    parts = []
    for im in imgs:
        m = re.match(r"\**(Fig\. [0-9a-z]+)", capt)
        alt = m.group(1) if m else "figure"
        parts.append(f"![{alt}]({im['src']})")
    if capt:
        capt = re.sub(r"^\*\*(Fig\. [^*]+)\*\*\s*", lambda m: "**" + m.group(1) + "** ", capt)
        parts.append(capt)
    return "\n\n".join(parts)


def walk(node, out):
    for c in node.children:
        if not isinstance(c, Tag):
            continue
        cls = c.get("class", [])
        if c.name == "section":
            out.append("\n---\n")
            walk(c, out)
        elif "rh" in cls or "folio" in cls or "frame" in cls or "bleed" in cls or "toc" in cls:
            continue
        elif "lbl" in cls:
            out.append(f"**{inline(c).strip()}**")
        elif "numeral" in cls:
            continue
        elif c.name == "h1":
            t = inline(c).strip().title() if inline(c).isupper() else inline(c).strip()
            out.append(f"# {t}")
        elif c.name == "h2":
            out.append(f"## {inline(c).strip()}")
        elif c.name == "h3":
            out.append(f"### {inline(c).strip()}")
        elif "eq" in cls:
            m = re.search(r"\\\[(.+?)\\\]", c.decode_contents(), re.S)
            if m:
                tex = html.unescape(m.group(1)).strip()
                out.append(f"$$\n{tex}\n$$")
        elif "q" in cls:
            out.append("> " + inline(c).strip())
        elif c.name == "p":
            t = inline(c).strip()
            if t:
                out.append(t)
        elif c.name in ("ul", "ol"):
            out.append("\n".join("- " + inline(li).strip() for li in c.find_all("li", recursive=False)))
        elif c.name == "table":
            out.append(table(c))
        elif c.name == "figure":
            out.append(figure(c))
        elif c.name == "figcaption":
            out.append(re.sub(r"^\*\*(Fig\. [^*]+)\*\*\s*", lambda m: "**" + m.group(1) + "** ", inline(c).strip()))
        elif c.name == "img":
            out.append(f"![]({c['src']})")
        else:
            walk(c, out)


def main():
    soup = BeautifulSoup((HERE / "report.html").read_text(encoding="utf-8"), "html.parser")
    out = []
    walk(soup.body, out)
    md = "\n\n".join(x for x in out if x.strip())
    md = re.sub(r"(\n\n---\n\n)+", "\n\n---\n\n", md)
    md = re.sub(r"\n{3,}", "\n\n", md)
    # cover and contents page become a title block
    head = ("# Re-entrant Auxetics\n\n*Cell models, simulation and a running midsole*\n\n"
            "MOS Project 2026–27 · Re-entrant strand\n\n"
            "**Contents:** [Summary](#summary) · [1 Re-entrant auxetics](#re-entrant-auxetics) · "
            "[2 Midsole optimisation](#midsole-optimisation) · [References](#references)\n\n")
    md = head + md[md.index("## Summary"):]
    (HERE / "report.md").write_text(md.strip() + "\n", encoding="utf-8")
    print("wrote report.md:", len(md.split()), "words")


if __name__ == "__main__":
    main()
