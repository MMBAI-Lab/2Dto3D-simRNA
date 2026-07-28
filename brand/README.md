# DansLab brand — vendored copy

Visual identity for DansLab (MMBAI), applied to everything this pipeline
publishes: the per-run and global reports under `results/` (`README_{EN,ES}.md`,
`REPORT_{EN,ES}.md`) and any future HTML/figure output.

**Upstream:** [MMBAI-Lab/branding-danslab](https://github.com/MMBAI-Lab/branding-danslab)
· vendored at commit `75b141f` (2026-07-22), same revision as the sibling
[2D-NAprediction](../../2D-NAprediction) repo.

The files here are a **copy**, not a submodule, so a clone of this repo builds
reports without needing access to the (private) brand repo. To refresh:

```bash
git clone --depth 1 git@github.com:MMBAI-Lab/branding-danslab.git /tmp/brand
cp /tmp/brand/{tokens.css,tokens.json,fonts.css} brand/
cp /tmp/brand/assets/danslab-logo-{dark,light}.svg brand/assets/
cp /tmp/brand/CLAUDE.md brand/BRAND-GUIDE.md
python -c "import cairosvg, io; from PIL import Image; \
Image.open(io.BytesIO(cairosvg.svg2png(url='brand/assets/danslab-logo-dark.svg', output_width=360)))\
.save('brand/assets/danslab-logo-report.png', dpi=(300,300))"
```

## Contents

| File | Use |
|---|---|
| `BRAND-GUIDE.md` | The rules (palette, type, layout, logo, voice). Read before restyling anything. |
| `tokens.json` | Canonical token values. Any generator that needs a color or font family reads it from here — never hardcode. |
| `tokens.css` | Same tokens as CSS custom properties, for HTML output. |
| `fonts.css` | Google Fonts `@import` for the IBM Plex family (needs network at view time). |
| `assets/danslab-logo-dark.svg` | Logo for light backgrounds (charcoal strokes, red atom). |
| `assets/danslab-logo-light.svg` | Logo for dark backgrounds. |
| `assets/danslab-logo-report.png` | 360 px @ 300 dpi raster of the dark logo — the report header mark. The dpi metadata is what makes pandoc place it at 1.2 in wide when a report is converted to DOCX; do not re-export without it. |

## How the reports use it

Every report under `results/` opens with the brand header block emitted by
[`scripts/report_brand.py`](../scripts/report_brand.py): the logo raster, the
`DANSLAB · MMBAI` kicker line, and the EN/ES language switch. Report generators
call `brand_header()` / `brand_footer()` rather than pasting the block, so a
change here propagates on the next build.

## Non-negotiables

- The logo's red atom `#E5322E` is fixed — never recolor the mark, rotate it,
  or add effects. Minimum width 88 px on screen / 18 mm in print.
- Red `#CE1B27` is the single accent: kickers, rules, numbers. Not for body text.
- Square corners, no drop shadows, no fonts or colors outside `tokens.json`.
- No emoji unless explicitly requested.
- Every public string ships EN + ES.
