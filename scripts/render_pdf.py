#!/usr/bin/env python3
"""Render a markdown file to PDF with local images embedded as base64 data URIs.

The make-pdf tool renders markdown through a headless Chromium, which refuses to load
local-file images (relative paths, absolute paths, and file:// URIs all come through blank —
only the alt-text/caption shows). This wrapper rewrites every *local* image reference into a
`data:` URI before rendering, so figures embed reliably. Remote (http/https) and existing
data: images are left untouched, and image paths resolve relative to the source markdown's
directory (so the committed .md keeps clean, readable relative paths).

Usage:
    scripts/render_pdf.py SRC.md OUT.pdf [extra make-pdf flags...]
    scripts/render_pdf.py paper.md paper.pdf --cover --toc --title "X" --author "Y"
"""
from __future__ import annotations

import base64
import os
import re
import subprocess
import sys
import tempfile

MIME = {
    ".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg",
    ".gif": "image/gif", ".svg": "image/svg+xml", ".webp": "image/webp",
}
IMG_RE = re.compile(r"(!\[[^\]]*\]\()([^)\s]+)(\s+\"[^\"]*\")?(\))")


def make_pdf_bin() -> str:
    cand = os.environ.get("MAKE_PDF_BIN") or os.path.expanduser(
        "~/.claude/skills/gstack/make-pdf/dist/pdf"
    )
    if not os.path.isfile(cand):
        sys.exit(f"make-pdf binary not found at {cand} (set MAKE_PDF_BIN)")
    return cand


def inline_images(md_text: str, base_dir: str) -> tuple[str, int]:
    n = 0

    def repl(m: re.Match) -> str:
        nonlocal n
        pre, src, title, close = m.group(1), m.group(2), m.group(3) or "", m.group(4)
        if src.startswith(("http://", "https://", "data:")):
            return m.group(0)
        path = src if os.path.isabs(src) else os.path.join(base_dir, src)
        if not os.path.isfile(path):
            print(f"  warn: image not found, left as-is: {src}", file=sys.stderr)
            return m.group(0)
        mime = MIME.get(os.path.splitext(path)[1].lower(), "image/png")
        b64 = base64.b64encode(open(path, "rb").read()).decode()
        n += 1
        # Emit raw HTML <img> (markdown passes it through) with max-width so wide matplotlib
        # figures scale to the page instead of overflowing the right margin. alt = caption.
        alt = re.match(r"!\[([^\]]*)\]", m.group(0)).group(1)
        return (f'<img src="data:{mime};base64,{b64}" alt="{alt}" '
                f'style="max-width:100%;height:auto;" />')

    return IMG_RE.sub(repl, md_text), n


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        sys.exit("usage: render_pdf.py SRC.md OUT.pdf [extra make-pdf flags...]")
    src, out, extra = argv[0], argv[1], argv[2:]
    base_dir = os.path.dirname(os.path.abspath(src))
    md_text = open(src).read()
    inlined, n = inline_images(md_text, base_dir)

    # write the temp md in the source dir so any other relative links still resolve
    fd, tmp = tempfile.mkstemp(suffix=".md", prefix=".render_", dir=base_dir)
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(inlined)
        cmd = [make_pdf_bin(), "generate", *extra, tmp, os.path.abspath(out)]
        proc = subprocess.run(cmd)
    finally:
        os.unlink(tmp)
    if proc.returncode == 0:
        print(f"Embedded {n} local image(s) -> {out}")
    return proc.returncode


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
