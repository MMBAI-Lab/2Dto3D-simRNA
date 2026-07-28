#!/usr/bin/env python3
"""
Compose the 2x2 summary PNG of the four centroids selected for MD: stitches
the already-ChimeraX-rendered centroid PNGs into a single figure with labels
(run + cluster rank + population) and a joint title.

Writes to results/APT-PF1/global_analysis/summary_md_selection.png.
"""
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

GLOBAL = Path("/home/pdans/Documents/2Dto3D-simRNA/results/APT-PF1/global_analysis")

SELECTION = [
    ("VFold2D",   "cutoff_12A/clust01", "DNA_as_RNA_VFold2D", 60.4),
    ("VFold2D",   "cutoff_12A/clust02", "DNA_as_RNA_VFold2D", 23.1),
    ("NUPACK4",   "cutoff_12A/clust01", "DNA_as_RNA_NUPACK4", 91.8),
    ("NA_as_DNA", "cutoff_12A/clust01", "NA_as_DNA_NUPACK4",  79.6),
]


def try_fonts():
    for cand in (
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "/usr/share/fonts/truetype/ubuntu/Ubuntu-B.ttf",
    ):
        if Path(cand).exists():
            return cand
    return None


def main():
    font_path = try_fonts()
    title_font = ImageFont.truetype(font_path, 32) if font_path else ImageFont.load_default()
    label_font = ImageFont.truetype(font_path, 22) if font_path else ImageFont.load_default()

    tiles = []
    for short, subpath, run_dir, pct in SELECTION:
        png = GLOBAL / run_dir / f"{subpath}.png"
        img = Image.open(png).convert("RGBA")
        tiles.append((short, subpath.split("/")[-1], pct, img))

    tile_w, tile_h = tiles[0][3].size
    pad_top = 80   # room for global title
    pad_label = 60 # room for per-panel label
    gap = 20

    canvas_w = 2 * tile_w + 3 * gap
    canvas_h = pad_top + 2 * (tile_h + pad_label + gap) + gap
    canvas = Image.new("RGBA", (canvas_w, canvas_h), (255, 255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    title = "APT-PF1 — 4 centroids selected for AMBER MD (cutoff 12 Å, > 15 % population)"
    draw.text((gap, 20), title, fill=(20, 20, 20), font=title_font)

    positions = [(0, 0), (1, 0), (0, 1), (1, 1)]
    for (col, row), (short, cluster, pct, img) in zip(positions, tiles):
        x = gap + col * (tile_w + gap)
        y = pad_top + row * (tile_h + pad_label + gap)
        canvas.paste(img, (x, y), img)
        label = f"{short} · {cluster} · {pct:.1f} %"
        draw.text((x, y + tile_h + 10), label, fill=(20, 20, 20), font=label_font)

    out = GLOBAL / "summary_md_selection.png"
    canvas.convert("RGB").save(out, "PNG")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
