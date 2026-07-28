# DansLab brand — instructions for Claude Code

When working in this project, apply the DansLab (MMBAI) visual identity. Use the
tokens in `tokens.css` / `tokens.json`; do not invent colors, fonts or spacing.

## Identity in one paragraph
DansLab is a molecular-modeling, bioinformatics & AI research group at Universidad
de la República (Salto, Uruguay). The visual system is editorial + computational:
IBM Plex type, a black/gray/white/red palette, a faint background grid, red corner
crosshairs on dark surfaces, and a molecule logo with one fixed red atom.

## Color
- Black `#000000` — dark surface (prefer full black, not off-black).
- Ink `#0B0B0C` — text on light.
- Graphite `#3A3B3F` · Ash `#9A9CA1` · Mist `#D9DADD` — neutral scale (body, captions, rules).
- Red `#CE1B27` — the single accent. Use sparingly (kickers, numbers, rules, active state).
- Dark red `#7D1015` — seals, hover, underlines.
- Paper `#F7F6F4` · White `#FFFFFF` — light surfaces.
- Logo red `#E5322E` — the mark's atom ONLY. Never use elsewhere and never change it.

## Type
- Headlines / titles: **IBM Plex Serif** (500–700).
- Body / running text: **IBM Plex Sans** (300–600).
- Labels, metadata, data, code: **IBM Plex Mono** (uppercase, letter-spaced for kickers).
- Load via `fonts.css` (Google Fonts).

## Layout system
- **Dual surface**: black for title & section-divider screens; white for content, data, documents.
- **Grid**: a 60px lattice underlies surfaces — line `rgba(255,255,255,.12)` on black, `rgba(0,0,0,.04)` on white.
- **Corner crosshairs**: red L-shaped marks in the corners of black surfaces.
- Corners are square — no border-radius on brand containers. No drop shadows on the logo.
- Kickers are mono, uppercase, letter-spacing ~.2em, in red.

## Logo
- Files: `assets/danslab-logo-dark.svg` (light bg) and `assets/danslab-logo-light.svg` (dark bg).
- The red atom is fixed. Do not recolor the mark, recolor the atom, rotate/skew, or add effects.
- Minimum width 88px on screen / 18mm in print. Clear space = the red atom's diameter on all sides.

## Voice
Precise, plain, unhyped. American English + Uruguayan/LATAM Spanish; both languages
always move together. Avoid marketing hype and dramatic one-liners.

## Don't
- No colors, fonts, gradients or rounded cards outside this system.
- No emoji unless explicitly requested.
- Don't restyle or animate the logo beyond the provided SVGs.
