"""misskappa hex logo (r-package/man/figures/logo.png), drawn on a 64x74 pixel
grid in the Sega Master System palette (each RGB channel one of 00/55/AA/FF),
upscaled with nearest neighbour. Run: python3 dev/logo/make_logo.py"""
import math
from pathlib import Path
from collections import deque
from PIL import Image

W, H, SCALE = 64, 74, 8
cx, cy = (W - 1) / 2, (H - 1) / 2
P = {  # SMS palette picks
    "K": (0, 0, 0), "w": (255, 255, 255),
    "g0": (0, 85, 0), "g1": (0, 170, 0), "g2": (85, 255, 85), "g3": (170, 255, 170),
    "p0": (170, 0, 85), "p1": (255, 85, 170), "p2": (255, 170, 255),
    "cy": (170, 255, 255), "or": (255, 170, 0),
}
img = {}

def inside_hex(x, y):
    dx, dy = abs(x - cx), abs(y - cy)
    return dx <= W / 2 - 0.5 and dy + dx / math.sqrt(3) <= H / 2 - 0.5

hexmask = {(x, y) for x in range(W) for y in range(H) if inside_hex(x, y)}
# depth from the hex edge (4-neighbour BFS) for a stepped pixel border
depth, q = {}, deque()
for (x, y) in hexmask:
    if any((x + a, y + b) not in hexmask for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1))):
        depth[(x, y)] = 0; q.append((x, y))
while q:
    x, y = q.popleft()
    for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        n = (x + a, y + b)
        if n in hexmask and n not in depth:
            depth[n] = depth[(x, y)] + 1; q.append(n)

bcol = ["K", "g2", "K"]  # stepped border: black, bright green, black
for p, d in depth.items():
    if d < len(bcol): img[p] = bcol[d]
    else:
        x, y = p
        # background: dark green, with a dithered lighter band toward the bottom
        img[p] = "g1" if (y > 59 and (x + y) % 2 == 0) or y > 64 else "g0"

# a few sparkles
for (x, y) in [(11, 25), (50, 18), (54, 38), (12, 56), (52, 57), (40, 9)]:
    if depth.get((x, y), 0) >= len(bcol) + 1:
        img[(x, y)] = "p2"
        for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            img[(x + a, y + b)] = "p1"

# --- creature (kappa yokai) -------------------------------------------------
spr = {}
fx, fy = 31.5, 34
def ell(x, y, ex, ey, rx, ry): return ((x - ex) / rx) ** 2 + ((y - ey) / ry) ** 2 <= 1
for x in range(W):
    for y in range(H):
        if ell(x, y, fx, fy, 11, 9.5): spr[(x, y)] = "g2"
        elif ell(x, y, fx, fy - 1.5, 13.5, 11.5) and y <= fy + 3: spr[(x, y)] = "g1"   # hair
for x in range(W):
    for y in range(H):
        if ell(x, y, fx, fy - 11, 7.5, 2.6):                                         # dish
            spr[(x, y)] = "w" if (y < fy - 11 and x < fx - 2) else "cy"
# hair fringe across the forehead
for x in range(19, 45):
    if ((x // 3) % 2 == 0 and spr.get((x, 25)) == "g2"): spr[(x, 25)] = "g1"
for x in range(18, 46):
    if spr.get((x, 24)) == "g2": spr[(x, 24)] = "g1"
for ex in (25, 37):                                                                   # eyes
    for y in range(31, 35):
        spr[(ex, y)] = spr[(ex + 1, y)] = "K"
    spr[(ex + 1, 31)] = "w"
for x in (22, 23, 40, 41): spr[(x, 36)] = "p1"                                        # cheeks
for x in range(29, 35): spr[(x, 37)] = "or"                                           # beak
for x in range(30, 34): spr[(x, 38)] = "or"

# outline
out = set()
for (x, y) in spr:
    for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        n = (x + a, y + b)
        if n not in spr: out.add(n)
for p in out:
    if p in hexmask and depth[p] >= len(bcol): img[p] = "K"
for p, c in spr.items():
    if p in hexmask: img[p] = c

# --- text: 10-row lowercase pixel font (0-1 ascender, 2-7 x-height, 8-9 descender)
F = {
 "m": [".....", ".....", "XXXX.", "X.X.X", "X.X.X", "X.X.X", "X.X.X", "X.X.X", ".....", "....."],
 "i": ["X", ".", "X", "X", "X", "X", "X", "X", ".", "."],
 "s": ["....", "....", ".XXX", "X...", "XX..", "..XX", "...X", "XXX.", "....", "...."],
 "k": ["X...", "X...", "X..X", "X.X.", "XX..", "XX..", "X.X.", "X..X", "....", "...."],
 "a": ["....", "....", ".XX.", "...X", ".XXX", "X..X", "X..X", ".XXX", "....", "...."],
 "p": ["....", "....", "XXX.", "X..X", "X..X", "X..X", "X..X", "XXX.", "X...", "X..."],
}
word = "misskappa"
tw = sum(len(F[c][0]) for c in word) + len(word) - 1
tx, ty = round(cx - tw / 2 + 0.5), 41
txt = {}
x = tx
for c in word:
    g = F[c]
    for r, row in enumerate(g):
        for j, ch in enumerate(row):
            if ch == "X": txt[(x + j, ty + r)] = "p2" if r <= 4 else "p1"
    x += len(g[0]) + 1
outl = set()
for (x, y) in txt:
    for a in (-1, 0, 1):
        for b in (-1, 0, 1, 2):          # 8-neighbour outline, one pixel deeper below
            n = (x + a, y + b)
            if n not in txt: outl.add(n)
for p in outl: img[p] = "K"
for p, c in txt.items(): img[p] = c

im = Image.new("RGBA", (W, H), (0, 0, 0, 0))
for (x, y), c in img.items():
    if (x, y) in hexmask: im.putpixel((x, y), P[c] + (255,))
out = Path(__file__).resolve().parents[2] / "r-package/man/figures/logo.png"
im.resize((W * SCALE, H * SCALE), Image.NEAREST).save(out)
