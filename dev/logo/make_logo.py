"""misskappa hex logo (r-package/man/figures/logo.png): a grim pixel kappa on a 64x74 grid, Sega Master System
palette (each RGB channel in {00,55,AA,FF}), nearest-neighbour upscale.
Run: python3 dev/logo/make_logo.py"""
import math
from collections import deque
from pathlib import Path

from PIL import Image

W, H, SCALE = 64, 74, 8
cx, cy = (W - 1) / 2, (H - 1) / 2
P = {
    "K": (0, 0, 0), "w": (255, 255, 255), "gy": (170, 170, 170), "dg": (85, 85, 85),
    "g0": (0, 85, 0), "g1": (0, 170, 0), "g2": (85, 255, 85),
    "sk": (85, 170, 85), "sd": (0, 85, 0), "sl": (170, 255, 85),     # skin, spots, skin highlight
    "bk": (170, 170, 85), "bl": (255, 255, 170),                      # beak / belly
    "h": (0, 0, 0), "hl": (85, 0, 85),                            # hair
    "wa": (85, 170, 255), "wd": (0, 85, 170),                         # dish water
    "ye": (255, 255, 85), "rd": (170, 0, 0), "mo": (85, 0, 0),        # eyes, mouth
    "sh": (170, 85, 0), "sx": (85, 0, 0),                             # shell
    "p0": (170, 0, 85), "p1": (255, 85, 170), "p2": (255, 170, 255),
}
img = {}
N4 = ((1, 0), (-1, 0), (0, 1), (0, -1))

def inside_hex(x, y):
    dx, dy = abs(x - cx), abs(y - cy)
    return dx <= W / 2 - 0.5 and dy + dx / math.sqrt(3) <= H / 2 - 0.5

hexmask = {(x, y) for x in range(W) for y in range(H) if inside_hex(x, y)}
depth, q = {}, deque()
for (x, y) in hexmask:
    if any((x + a, y + b) not in hexmask for a, b in N4):
        depth[(x, y)] = 0; q.append((x, y))
while q:
    x, y = q.popleft()
    for a, b in N4:
        n = (x + a, y + b)
        if n in hexmask and n not in depth:
            depth[n] = depth[(x, y)] + 1; q.append(n)
bcol = ["K", "g2", "K"]
for p, d in depth.items():
    img[p] = bcol[d] if d < len(bcol) else "g0"

# reeds along both sides
for x0, h0, lean in [(6, 30, 0), (9, 38, 1), (12, 26, 0), (51, 34, 0), (54, 40, -1), (57, 28, 0), (15, 20, 1), (48, 22, -1)]:
    for k in range(h0):
        x, y = x0 + (lean * k) // 9, 70 - k
        if (x, y) in depth and depth[(x, y)] >= 3: img[(x, y)] = "g1"

# --- the kappa -------------------------------------------------------------
spr = {}
fx = 31.5
def ell(x, y, ex, ey, rx, ry): return ((x - ex) / rx) ** 2 + ((y - ey) / ry) ** 2 <= 1
def tri(x, y, a, b, c):
    def s(p, q, r): return (p[0] - r[0]) * (q[1] - r[1]) - (q[0] - r[0]) * (p[1] - r[1])
    d1, d2, d3 = s((x, y), a, b), s((x, y), b, c), s((x, y), c, a)
    return not ((d1 < 0 or d2 < 0 or d3 < 0) and (d1 > 0 or d2 > 0 or d3 > 0))
G = [(x, y) for x in range(W) for y in range(H)]

# body: hunched shoulders + shell rim, behind the title
for (x, y) in G:
    if ell(x, y, fx, 64, 30, 22) and y >= 42:
        spr[(x, y)] = "sh" if (x + 2 * y) % 7 else "sx"
    if ell(x, y, fx, 64, 25, 20) and y >= 44:
        spr[(x, y)] = "sk"
        if abs(x - fx) <= 6 and y >= 50: spr[(x, y)] = "bl" if y % 3 else "bk"   # belly plates
# neck
for (x, y) in G:
    if 38 <= y <= 46 and abs(x - fx) <= 6: spr[(x, y)] = "sk"
# ears: long, pointed, swept up and out
for (x, y) in G:
    if tri(x, y, (fx - 9, 25), (fx - 9, 32), (fx - 19, 19)) or tri(x, y, (fx + 9, 25), (fx + 9, 32), (fx + 19, 19)):
        spr[(x, y)] = "sk"
# head
for (x, y) in G:
    if ell(x, y, fx, 28, 11.5, 12.5): spr[(x, y)] = "sk"
# skin spots
for p in [(25, 21), (38, 22), (23, 31), (40, 30), (27, 34), (36, 35), (29, 43), (35, 41), (19, 50), (44, 49), (14, 55), (49, 56)]:
    if spr.get(p) == "sk": spr[p] = "sd"
# bowl-cut hair with ragged fringe and side locks
for (x, y) in G:
    if ell(x, y, fx, 21, 12.5, 10.5) and (y <= 23 or (y <= 30 and abs(x - fx) >= 9)):
        spr[(x, y)] = "h"
for x in range(int(fx) - 10, int(fx) + 12):
    if x % 3 != 1: spr[(x, 24)] = "h"
for (x, y) in G:   # hair highlight streaks
    if (x, y) in {(23, 17), (24, 16), (25, 15), (22, 19), (23, 18), (37, 15), (38, 16), (39, 17)}: spr[(x, y)] = "hl"
# dish on the crown
for (x, y) in G:
    if ell(x, y, fx, 12, 6.5, 2.4): spr[(x, y)] = "gy"
    if ell(x, y, fx, 12, 4.6, 1.3): spr[(x, y)] = "wa"
spr[(29, 12)] = "w"; spr[(30, 12)] = "w"; spr[(34, 13)] = "wd"
# eyes: yellow, glaring sideways, under a heavy slanted brow
for ex, look in ((25, 1), (35, 1)):
    for x in range(ex, ex + 4):
        for y in (26, 27, 28): spr[(x, y)] = "ye"
    spr[(ex + 1 + look, 27)] = "K"; spr[(ex + 1 + look, 28)] = "K"; spr[(ex + 1 + look - 1, 27)] = "rd"
inner_l, inner_r = 28, 35
for x in range(24, 30): spr[(x, 25 + (x - 24) // 3)] = "K"          # brows slope down to the middle
for x in range(34, 40): spr[(x, 25 + (39 - x) // 3)] = "K"
# beak: pale, hooked, gaping with teeth
for (x, y) in G:
    BW = {29: 1, 30: 1.5, 31: 2.5, 32: 3.5, 33: 4.5, 34: 5, 35: 5, 36: 5, 37: 5, 38: 4, 39: 3, 40: 1.5}
    if y in BW and abs(x - fx) <= BW[y]: spr[(x, y)] = "bk"
for x in (30, 33): spr[(x, 32)] = "K"                               # nostrils
for (x, y) in G:
    if 34 <= y <= 37 and abs(x - fx) <= 4: spr[(x, y)] = "mo"
for x in range(28, 36):                                            # top and bottom teeth
    spr[(x, 34)] = "w" if x % 2 == 0 else "mo"
    spr[(x, 37)] = "w" if x % 2 == 1 else "mo"
# outline
for (x, y) in list(spr):
    for a, b in N4:
        n = (x + a, y + b)
        if n not in spr and n in depth and depth[n] >= 3: img[n] = "K"
for p, c in spr.items():
    if p in depth and depth[p] >= 3: img[p] = c

# --- title -----------------------------------------------------------------
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
tx, ty = round(cx - tw / 2 + 0.5), 47
txt, x = {}, tx
for c in word:
    for r, row in enumerate(F[c]):
        for j, ch in enumerate(row):
            if ch == "X": txt[(x + j, ty + r)] = "p2" if r <= 4 else "p1"
    x += len(F[c][0]) + 1
for (x, y) in list(txt):
    for a in (-1, 0, 1):
        for b in (-1, 0, 1, 2):
            n = (x + a, y + b)
            if n not in txt: img[n] = "K"
img.update(txt)

im = Image.new("RGBA", (W, H), (0, 0, 0, 0))
for p, c in img.items():
    if p in hexmask: im.putpixel(p, P[c] + (255,))
out = Path(__file__).resolve().parents[2] / "r-package/man/figures/logo.png"
im.resize((W * SCALE, H * SCALE), Image.NEAREST).save(out)
