#!/usr/bin/env python3
"""Build branded DOCX reports from the Markdown reports under results/.

Every report ships as a bilingual pair. This repo has two basenames:
`README_{EN,ES}.md` for per-run substate summaries and `REPORT_{EN,ES}.md`
for cross-run analyses (see CLAUDE.md). Both are converted here into a DOCX
styled with the DansLab (MMBAI) brand system: IBM Plex type, the
black/graphite/red palette, square corners, red section rules, brand footer.

The styling lives in a pandoc reference document generated on the fly from
`brand/tokens.json` — there is no binary template committed to the repo, so a
change to the brand tokens propagates on the next build.

Ported from the sibling 2D-NAprediction repo; keep the two in sync if the
brand furniture changes.

Usage:
    python3 scripts/build_report.py                       # every report pair under results/
    python3 scripts/build_report.py results/comercialApt  # one set (all reports, both langs)
    python3 scripts/build_report.py results/comercialApt/REPORT_ES.md
    python3 scripts/build_report.py --reference-out /tmp/danslab-reference.docx
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from report_brand import BRAND_TOKENS, REPO_ROOT  # noqa: E402

RESULTS = REPO_ROOT / "results"

#: Report basenames, in the order they should be built. Anything matching
#: <BASE>_{EN,ES}.md under results/ is a report; nothing else is.
REPORT_BASES = ("REPORT", "README")

W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
FOOTER_REL_ID = "rId90"
FOOTER_TEXT = "DansLab · MMBAI — Universidad de la República, Salto, Uruguay"

# A4 portrait, 20 mm vertical / 16 mm horizontal margins, in twips.
# The narrow side margins are what let a 74 nt dot-bracket stay on one line.
PAGE_W, PAGE_H = 11906, 16838
MARGIN_V, MARGIN_H = 1134, 907
TEXT_WIDTH = PAGE_W - 2 * MARGIN_H


def _tokens() -> dict:
    return json.loads(BRAND_TOKENS.read_text())


class Brand:
    """Brand tokens in the units Word wants (hex without '#', half-points)."""

    def __init__(self) -> None:
        color = _tokens()["color"]
        self.black = color["black"].lstrip("#").upper()
        self.ink = color["ink"].lstrip("#").upper()
        self.graphite = color["graphite"].lstrip("#").upper()
        self.slate = color["slate"].lstrip("#").upper()
        self.ash = color["ash"].lstrip("#").upper()
        self.mist = color["mist"].lstrip("#").upper()
        self.paper = color["paper"].lstrip("#").upper()
        self.red = color["red"].lstrip("#").upper()
        self.red_dark = color["redDark"].lstrip("#").upper()
        # tokens.json ships CSS font stacks; the DOCX needs the family name.
        self.serif = _first_family(_tokens()["font"]["serif"])
        self.sans = _first_family(_tokens()["font"]["sans"])
        self.mono = _first_family(_tokens()["font"]["mono"])


def _first_family(stack: str) -> str:
    return stack.split(",")[0].strip().strip('"')


def _fonts(family: str) -> str:
    return (f'<w:rFonts w:ascii="{family}" w:hAnsi="{family}" '
            f'w:eastAsia="{family}" w:cs="{family}"/>')


def _border(edge: str, size: int, color: str, space: int = 4) -> str:
    return f'<w:{edge} w:val="single" w:sz="{size}" w:space="{space}" w:color="{color}"/>'


def _style(style_id: str, name: str, *, kind: str = "paragraph",
           based_on: str | None = "Normal", next_style: str | None = None,
           ppr: str = "", rpr: str = "", default: bool = False,
           custom: bool = False) -> str:
    parts = [f'<w:style w:type="{kind}"'
             + (' w:default="1"' if default else "")
             + f' w:customStyle="{1 if custom else 0}" w:styleId="{style_id}">',
             f'<w:name w:val="{name}"/>']
    if based_on:
        parts.append(f'<w:basedOn w:val="{based_on}"/>')
    if next_style:
        parts.append(f'<w:next w:val="{next_style}"/>')
    parts.append("<w:qFormat/>")
    if ppr:
        parts.append(f"<w:pPr>{ppr}</w:pPr>")
    if rpr:
        parts.append(f"<w:rPr>{rpr}</w:rPr>")
    parts.append("</w:style>")
    return "".join(parts)


def _heading(level: int, family: str, size: int, color: str, *,
             before: int, after: int, bottom_rule: str = "") -> str:
    ppr = "<w:keepNext/>"
    if bottom_rule:
        ppr += f"<w:pBdr>{bottom_rule}</w:pBdr>"
    ppr += f'<w:spacing w:before="{before}" w:after="{after}"/>'
    rpr = _fonts(family) + f'<w:b/><w:color w:val="{color}"/>' \
        + f'<w:sz w:val="{size}"/><w:szCs w:val="{size}"/>'
    return _style(f"Heading{level}", f"heading {level}", next_style="BodyText",
                  ppr=ppr, rpr=rpr)


def brand_styles(b: Brand) -> dict[str, str]:
    """Full replacements for the pandoc styles we restyle, keyed by styleId."""
    sans, serif, mono = b.sans, b.serif, b.mono
    styles = {
        "Normal": _style(
            "Normal", "Normal", based_on=None, default=True,
            ppr='<w:spacing w:before="0" w:after="140" w:line="264" w:lineRule="auto"/>',
            rpr=_fonts(sans) + f'<w:color w:val="{b.ink}"/><w:sz w:val="20"/><w:szCs w:val="20"/>'),

        # H1 carries the report title: serif, ink, red rule underneath.
        "Heading1": _heading(1, serif, 34, b.ink, before=0, after=200,
                             bottom_rule=_border("bottom", 16, b.red, space=6)),
        "Heading2": _heading(2, serif, 26, b.ink, before=360, after=140,
                             bottom_rule=_border("bottom", 6, b.mist, space=4)),
        "Heading3": _heading(3, serif, 22, b.graphite, before=280, after=100),
        "Heading4": _heading(4, sans, 20, b.graphite, before=240, after=80),
        "Heading5": _heading(5, sans, 19, b.slate, before=200, after=80),
        "Heading6": _heading(6, sans, 19, b.slate, before=200, after=80),

        "Title": _style(
            "Title", "Title", next_style="Subtitle",
            ppr='<w:pBdr>' + _border("bottom", 16, b.red, space=6) + '</w:pBdr>'
                '<w:spacing w:before="0" w:after="160"/>',
            rpr=_fonts(serif) + f'<w:b/><w:color w:val="{b.ink}"/><w:sz w:val="48"/><w:szCs w:val="48"/>'),
        "Subtitle": _style(
            "Subtitle", "Subtitle", next_style="BodyText",
            ppr='<w:spacing w:before="0" w:after="240"/>',
            rpr=_fonts(sans) + f'<w:color w:val="{b.slate}"/><w:sz w:val="22"/><w:szCs w:val="22"/>'),
        "Author": _style(
            "Author", "Author", next_style="BodyText",
            rpr=_fonts(mono) + f'<w:color w:val="{b.slate}"/><w:sz w:val="17"/><w:szCs w:val="17"/>'),
        "Date": _style(
            "Date", "Date", next_style="BodyText",
            rpr=_fonts(mono) + f'<w:color w:val="{b.slate}"/><w:sz w:val="17"/><w:szCs w:val="17"/>'),

        # Code: dot-brackets must stay on one line, hence the small mono size.
        "VerbatimChar": _style(
            "VerbatimChar", "Verbatim Char", kind="character",
            based_on="DefaultParagraphFont",
            rpr=_fonts(mono) + f'<w:color w:val="{b.ink}"/><w:sz w:val="16"/><w:szCs w:val="16"/>'),
        "SourceCode": _style(
            "SourceCode", "Source Code", next_style="BodyText",
            ppr='<w:pBdr>' + _border("left", 12, b.mist, space=6) + '</w:pBdr>'
                f'<w:shd w:val="clear" w:color="auto" w:fill="{b.paper}"/>'
                '<w:spacing w:before="140" w:after="140" w:line="240" w:lineRule="auto"/>'
                '<w:ind w:left="120"/><w:contextualSpacing w:val="0"/>',
            rpr=_fonts(mono) + f'<w:color w:val="{b.ink}"/><w:sz w:val="16"/><w:szCs w:val="16"/>'),

        # Blockquotes are used for callouts (confidentiality, tool caveats).
        "BlockText": _style(
            "BlockText", "Block Text",
            ppr='<w:pBdr>' + _border("left", 12, b.red, space=8) + '</w:pBdr>'
                '<w:ind w:left="200"/><w:spacing w:before="160" w:after="160"/>',
            rpr=_fonts(sans) + f'<w:color w:val="{b.graphite}"/><w:sz w:val="19"/><w:szCs w:val="19"/>'),

        "Caption": _style(
            "Caption", "Caption",
            ppr='<w:spacing w:before="80" w:after="240"/>',
            rpr=_fonts(sans) + f'<w:color w:val="{b.slate}"/><w:sz w:val="17"/><w:szCs w:val="17"/>'),
        "ImageCaption": _style(
            "ImageCaption", "Image Caption", based_on="Caption"),
        "TableCaption": _style(
            "TableCaption", "Table Caption", based_on="Caption"),
        "Figure": _style(
            "Figure", "Figure", next_style="ImageCaption",
            ppr='<w:spacing w:before="200" w:after="80"/>'),
        "CaptionedFigure": _style(
            "CaptionedFigure", "Captioned Figure", based_on="Figure"),

        "FootnoteText": _style(
            "FootnoteText", "footnote text",
            rpr=_fonts(sans) + f'<w:color w:val="{b.slate}"/><w:sz w:val="16"/><w:szCs w:val="16"/>'),
        "Hyperlink": _style(
            "Hyperlink", "Hyperlink", kind="character",
            based_on="DefaultParagraphFont",
            rpr=f'<w:color w:val="{b.red_dark}"/><w:u w:val="single"/>'),

        # Editorial table: horizontal rules only.
        "Table": _style(
            "Table", "Table", kind="table", based_on="TableNormal", default=True,
            ppr="", rpr=_fonts(sans) + f'<w:color w:val="{b.ink}"/><w:sz w:val="17"/><w:szCs w:val="17"/>'),
    }
    # The table style needs <w:tblPr>, which _style() does not model.
    styles["Table"] = styles["Table"].replace(
        "</w:style>",
        "<w:tblPr><w:tblBorders>"
        + _border("top", 8, b.ink, space=0)
        + _border("bottom", 8, b.ink, space=0)
        + _border("insideH", 4, b.mist, space=0)
        + "</w:tblBorders>"
        "<w:tblCellMar>"
        '<w:top w:w="60" w:type="dxa"/><w:left w:w="80" w:type="dxa"/>'
        '<w:bottom w:w="60" w:type="dxa"/><w:right w:w="80" w:type="dxa"/>'
        "</w:tblCellMar></w:tblPr></w:style>")
    return styles


def extra_styles(b: Brand) -> str:
    """Custom styles the Lua filter attaches to the report boilerplate."""
    kicker = _style(
        "Kicker", "Kicker", kind="character", based_on="DefaultParagraphFont",
        custom=True,
        rpr=_fonts(b.mono) + f'<w:color w:val="{b.red}"/><w:caps/>'
            '<w:spacing w:val="30"/><w:sz w:val="15"/><w:szCs w:val="15"/>')
    meta = _style(
        "Meta", "Meta", kind="character", based_on="DefaultParagraphFont",
        custom=True,
        rpr=_fonts(b.mono) + f'<w:color w:val="{b.slate}"/>'
            '<w:sz w:val="16"/><w:szCs w:val="16"/>')
    # Dot-brackets run to 70+ characters and must not wrap inside a table
    # cell, so they get a size of their own (6.5 pt) below inline code (8 pt).
    dotbracket = _style(
        "DotBracket", "Dot Bracket", kind="character",
        based_on="DefaultParagraphFont", custom=True,
        rpr=_fonts(b.mono) + f'<w:color w:val="{b.ink}"/>'
            '<w:sz w:val="13"/><w:szCs w:val="13"/>')
    footer = _style(
        "Footer", "footer",
        ppr='<w:pBdr>' + _border("top", 4, b.mist, space=6) + '</w:pBdr>'
            f'<w:tabs><w:tab w:val="right" w:pos="{TEXT_WIDTH}"/></w:tabs>'
            '<w:spacing w:before="0" w:after="0"/>',
        rpr=_fonts(b.mono) + f'<w:color w:val="{b.slate}"/>'
            '<w:sz w:val="14"/><w:szCs w:val="14"/>')
    return kicker + meta + dotbracket + footer


def footer_xml(b: Brand) -> str:
    run_props = (f"<w:rPr>{_fonts(b.mono)}"
                 f'<w:color w:val="{b.slate}"/><w:sz w:val="14"/><w:szCs w:val="14"/></w:rPr>')
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        f'<w:ftr xmlns:w="{W_NS}" xmlns:r="{R_NS}">'
        '<w:p><w:pPr><w:pStyle w:val="Footer"/></w:pPr>'
        f'<w:r>{run_props}<w:t xml:space="preserve">{FOOTER_TEXT}</w:t></w:r>'
        f'<w:r>{run_props}<w:tab/></w:r>'
        f'<w:r>{run_props}<w:fldChar w:fldCharType="begin"/></w:r>'
        f'<w:r>{run_props}<w:instrText xml:space="preserve"> PAGE </w:instrText></w:r>'
        f'<w:r>{run_props}<w:fldChar w:fldCharType="separate"/></w:r>'
        f'<w:r>{run_props}<w:t>1</w:t></w:r>'
        f'<w:r>{run_props}<w:fldChar w:fldCharType="end"/></w:r>'
        '</w:p></w:ftr>')


def sect_pr() -> str:
    return (
        f'<w:sectPr><w:footerReference w:type="default" r:id="{FOOTER_REL_ID}"/>'
        f'<w:pgSz w:w="{PAGE_W}" w:h="{PAGE_H}"/>'
        f'<w:pgMar w:top="{MARGIN_V}" w:right="{MARGIN_H}" w:bottom="{MARGIN_V}" '
        f'w:left="{MARGIN_H}" w:header="708" w:footer="454" w:gutter="0"/>'
        '<w:cols w:space="708"/><w:docGrid w:linePitch="360"/></w:sectPr>')


def patch_sect_pr(xml: str) -> str:
    """Swap the reference document's page setup for ours.

    pandoc reuses the *first* sectPr it finds in the reference document, so
    the empty `<w:sectPr />` that ships with the default reference has to be
    replaced rather than appended to.
    """
    import re

    pattern = re.compile(r"<w:sectPr\s*/>|<w:sectPr(?: [^>]*)?>.*?</w:sectPr>", re.S)
    if pattern.search(xml):
        return pattern.sub(lambda _m: sect_pr(), xml, count=1)
    return xml.replace("</w:body>", sect_pr() + "</w:body>")


def patch_styles(xml: str, b: Brand) -> str:
    import re

    for style_id, replacement in brand_styles(b).items():
        pattern = re.compile(
            r'<w:style [^>]*w:styleId="' + re.escape(style_id) + r'"[^>]*>.*?</w:style>',
            re.S)
        xml, n = pattern.subn(lambda _m, r=replacement: r, xml, count=1)
        if not n:
            xml = xml.replace("</w:styles>", replacement + "</w:styles>")
    xml = xml.replace("</w:styles>", extra_styles(b) + "</w:styles>")
    # Anything still resolving through the theme lands on Plex too.
    return xml


def patch_theme(xml: str, b: Brand) -> str:
    import re

    xml = re.sub(r'(<a:majorFont><a:latin typeface=")[^"]*(")',
                 lambda m: m.group(1) + b.serif + m.group(2), xml, count=1)
    xml = re.sub(r'(<a:minorFont><a:latin typeface=")[^"]*(")',
                 lambda m: m.group(1) + b.sans + m.group(2), xml, count=1)
    xml = re.sub(r'<a:accent1><a:srgbClr val="[^"]*" ?/></a:accent1>',
                 f'<a:accent1><a:srgbClr val="{b.red}"/></a:accent1>', xml, count=1)
    return xml


def build_reference_docx(dest: Path) -> Path:
    """Generate the DansLab pandoc reference document at `dest`."""
    b = Brand()
    base = subprocess.run(["pandoc", "--print-default-data-file", "reference.docx"],
                          check=True, capture_output=True).stdout
    with tempfile.TemporaryDirectory() as td:
        base_path = Path(td) / "base.docx"
        base_path.write_bytes(base)
        dest.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(base_path) as zin, \
                zipfile.ZipFile(dest, "w", zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                data = zin.read(item.filename)
                if item.filename == "word/styles.xml":
                    data = patch_styles(data.decode("utf-8"), b).encode("utf-8")
                elif item.filename == "word/theme/theme1.xml":
                    data = patch_theme(data.decode("utf-8"), b).encode("utf-8")
                elif item.filename == "word/document.xml":
                    data = patch_sect_pr(data.decode("utf-8")).encode("utf-8")
                elif item.filename == "word/_rels/document.xml.rels":
                    rel = (f'<Relationship Id="{FOOTER_REL_ID}" '
                           f'Type="{R_NS}/footer" Target="footer1.xml"/>')
                    data = data.decode("utf-8").replace(
                        "</Relationships>", rel + "</Relationships>").encode("utf-8")
                elif item.filename == "[Content_Types].xml":
                    override = ('<Override PartName="/word/footer1.xml" ContentType='
                                '"application/vnd.openxmlformats-officedocument.'
                                'wordprocessingml.footer+xml"/>')
                    data = data.decode("utf-8").replace(
                        "</Types>", override + "</Types>").encode("utf-8")
                zout.writestr(item, data)
            zout.writestr("word/footer1.xml", footer_xml(b))
    return dest


def build_docx(md: Path, reference: Path, lua_filter: Path) -> Path:
    out = md.with_suffix(".docx")
    cmd = [
        "pandoc", str(md), "-o", str(out),
        "--from", "markdown", "--to", "docx",
        "--reference-doc", str(reference),
        "--lua-filter", str(lua_filter),
        "--resource-path", f"{md.parent}:{REPO_ROOT}",
        "--standalone",
    ]
    subprocess.run(cmd, check=True)
    return out


def collect_targets(paths: list[Path]) -> list[Path]:
    if not paths:
        paths = [RESULTS]
    targets: list[Path] = []
    for p in paths:
        if p.is_file():
            targets.append(p)
        elif p.is_dir():
            for base in REPORT_BASES:
                for lang in ("EN", "ES"):
                    targets.extend(sorted(p.glob(f"**/{base}_{lang}.md")))
        else:
            raise SystemExit(f"not found: {p}")
    return targets


def warn_unpaired(targets: list[Path]) -> int:
    """Reports are bilingual by policy — flag any half-translated pair."""
    missing = 0
    for md in targets:
        base = md.name.rsplit("_", 1)[0]
        other = "ES" if md.name.endswith("_EN.md") else "EN"
        sibling = md.parent / f"{base}_{other}.md"
        if not sibling.exists():
            missing += 1
            print(f"[warn] {md} has no {sibling.name} — reports must ship EN + ES")
    return missing


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="*", type=Path,
                    help="REPORT_*.md / README_*.md files, or directories to scan "
                         "(default: results/)")
    ap.add_argument("--reference-out", type=Path,
                    help="Only write the branded pandoc reference document here "
                         "(for inspection) and exit.")
    args = ap.parse_args()

    if shutil.which("pandoc") is None:
        raise SystemExit("pandoc not found — install it (apt install pandoc)")

    if args.reference_out:
        print(f"Wrote {build_reference_docx(args.reference_out)}")
        return 0

    targets = collect_targets(args.paths)
    if not targets:
        print("No report pairs found (looked for "
              + " / ".join(f"{b}_{{EN,ES}}.md" for b in REPORT_BASES) + ").")
        return 0
    warn_unpaired(targets)

    lua_filter = Path(__file__).resolve().parent / "danslab_docx.lua"
    with tempfile.TemporaryDirectory() as td:
        reference = build_reference_docx(Path(td) / "danslab-reference.docx")
        for md in targets:
            out = build_docx(md, reference, lua_filter)
            print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
