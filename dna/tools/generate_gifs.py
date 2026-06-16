#!/usr/bin/env python3
from __future__ import annotations

import math
import shutil
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "assets"
TMP = ROOT / "tmp" / "gif-frames"
W, H = 960, 540
BG = "#050608"


def svg_text(x: float, y: float, text: str, color: str, size: int = 22, weight: int = 700) -> str:
    return (
        f'<text x="{x:.1f}" y="{y:.1f}" fill="{color}" font-family="monospace" '
        f'font-size="{size}" font-weight="{weight}" text-anchor="middle" '
        f'dominant-baseline="middle">{text}</text>'
    )


def frame(path: Path, body: str) -> None:
    path.write_text(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">'
        f'<rect width="100%" height="100%" fill="{BG}"/>{body}</svg>',
        encoding="utf-8",
    )


def render_helix() -> None:
    frames = TMP / "helix"
    frames.mkdir(parents=True, exist_ok=True)
    seq = "ATGGCTGACGAATTCGGTCACTGCAAGTTGACCGTTGGTAC"
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    colors = {"A": "#ff5a5f", "C": "#5b8cff", "G": "#ffd166", "T": "#5ee37d"}
    rows = 24
    for f in range(36):
        theta = f * 0.22
        parts = [svg_text(W / 2, 34, "DNA codon helix", "#8b949e", 20, 500)]
        for y in range(rows):
            t = y * 0.42 + theta
            z = math.sin(t)
            x1 = W / 2 + math.cos(t) * 205
            x2 = W / 2 - math.cos(t) * 205
            yy = 72 + y * 18
            base = seq[(f // 2 + y) % len(seq)]
            pair = pairs[base]
            if z >= 0:
                front_x, back_x, front, back = x1, x2, base, pair
            else:
                front_x, back_x, front, back = x2, x1, pair, base
            parts.append(
                f'<line x1="{back_x:.1f}" y1="{yy:.1f}" x2="{front_x:.1f}" y2="{yy:.1f}" '
                f'stroke="#30363d" stroke-width="2"/>'
            )
            parts.append(svg_text(back_x, yy, back, colors[back], 24, 700))
            parts.append(svg_text(front_x, yy, front, colors[front], 24, 900))
        frame(frames / f"{f:03}.svg", "".join(parts))
    make_gif(frames, ASSETS / "helix.gif", delay=6)


def render_matrix() -> None:
    frames = TMP / "matrix"
    frames.mkdir(parents=True, exist_ok=True)
    bases = "ATCG"
    colors = {"A": "#ff5a5f", "C": "#5b8cff", "G": "#ffd166", "T": "#5ee37d"}
    cols = 58
    for f in range(36):
        parts = [svg_text(W / 2, 34, "ATCG matrix rain", "#8b949e", 20, 500)]
        for col in range(cols):
            x = 28 + col * 16
            head = (f * (1 + col % 3) + col * 11) % 34
            length = 7 + (col * 5) % 12
            for i in range(length):
                y = 68 + (head - i) * 15
                if y < 58 or y > H - 24:
                    continue
                base = bases[(col * 7 + i * 3 + f) % 4]
                opacity = 1.0 if i == 0 else max(0.25, 0.9 - i * 0.07)
                parts.append(
                    f'<text x="{x:.1f}" y="{y:.1f}" fill="{colors[base]}" fill-opacity="{opacity:.2f}" '
                    f'font-family="monospace" font-size="22" font-weight="800" text-anchor="middle">{base}</text>'
                )
        frame(frames / f"{f:03}.svg", "".join(parts))
    make_gif(frames, ASSETS / "matrix.gif", delay=5)


def make_gif(frames: Path, out: Path, delay: int) -> None:
    svgs = sorted(str(path) for path in frames.glob("*.svg"))
    subprocess.run(
        ["magick", "-delay", str(delay), "-loop", "0", *svgs, str(out)],
        check=True,
        cwd=ROOT,
    )


def main() -> None:
    if shutil.which("magick") is None:
        raise SystemExit("missing ImageMagick 'magick'")
    if TMP.exists():
        shutil.rmtree(TMP)
    ASSETS.mkdir(exist_ok=True)
    render_helix()
    render_matrix()


if __name__ == "__main__":
    main()
