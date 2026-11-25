import csv, json, sys, hashlib
from PIL import Image, ImageDraw, ImageFont
proj, exports, out = sys.argv[1], sys.argv[2], sys.argv[3]
img = Image.open(f"{proj}/Screenshot_1.png").convert("RGB")
d = json.load(open(f"{proj}/wpd.json"))
pts = [p for s in d["datasetColl"] if s["name"] == "exp_small" for p in s["data"]]
rows = list(csv.reader(open(f"{exports}/exp_small.csv")))[1:]
assert len(pts) == 9 and len(rows) == 9
for p, r in zip(pts, rows):  # stored order == export order
    assert abs(p["x"] - float(r[2])) < 1e-9 and abs(p["y"] - float(r[3])) < 1e-9
S = 2
big = img.resize((img.width * S, img.height * S), Image.LANCZOS)
dr = ImageDraw.Draw(big)
try:
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 22)
    fsm = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 16)
except OSError:
    font = fsm = ImageFont.load_default()
COL = (230, 0, 160)
# label offsets (in 1x px) chosen to avoid crowding at low k; leader line drawn to ring edge
off = {1: (-45, -40), 2: (-30, -62), 3: (30, 34), 4: (10, 38), 5: (10, 38),
       6: (10, 38), 7: (-30, -36), 8: (12, 34), 9: (-34, -30)}
R = 11  # ring radius, 1x px (marker radius is about 6 px): ring sits outside the marker
for i, p in enumerate(pts, 1):
    x, y = p["x"] * S, p["y"] * S
    dr.ellipse([x - R * S, y - R * S, x + R * S, y + R * S], outline=COL, width=2)
    ox, oy = off[i]
    lx, ly = x + ox * S, y + oy * S
    import math
    dx, dy = lx - x, ly - y; L = math.hypot(dx, dy)
    dr.line([x + dx / L * R * S, y + dy / L * R * S, lx, ly], fill=COL, width=2)
    t = str(i); bb = dr.textbbox((0, 0), t, font=font)
    w, h = bb[2] - bb[0], bb[3] - bb[1]
    dr.rectangle([lx - w / 2 - 5, ly - h / 2 - 6, lx + w / 2 + 5, ly + h / 2 + 6], fill="white", outline=COL, width=2)
    dr.text((lx - w / 2, ly - h / 2 - bb[1]), t, fill=COL, font=font)
canvas = Image.new("RGB", (big.width, big.height + 34), "white"); canvas.paste(big, (0, 0)); big = canvas; dr = ImageDraw.Draw(big)
dr.text((10, big.height - 26), "S3.5 visual check: exp_small points 1-9 in click order (rings drawn outside markers; coordinates from wpd.json)", fill=COL, font=fsm)
big.save(f"{out}/s35_overview_numbered.png")
# zoom sheet: raw pixels, 6x nearest, ticks with a gap so the marker itself is untouched
Z, H = 6, 22
tiles = []
for i, p in enumerate(pts, 1):
    cx, cy = round(p["x"]), round(p["y"])
    crop = img.crop((cx - H, cy - H, cx + H, cy + H)).resize((2 * H * Z, 2 * H * Z), Image.NEAREST)
    t = Image.new("RGB", (crop.width, crop.height + 30), "white"); t.paste(crop, (0, 30))
    td = ImageDraw.Draw(t)
    mx, my = (p["x"] - (cx - H)) * Z, (p["y"] - (cy - H)) * Z + 30
    g, l = 12 * Z, 6 * Z
    for a, b, c, e in [(mx - g - l, my, mx - g, my), (mx + g, my, mx + g + l, my),
                       (mx, my - g - l, mx, my - g), (mx, my + g, mx, my + g + l)]:
        td.line([a, b, c, e], fill=COL, width=3)
    td.text((6, 4), f"{i}: px={p['x']:.1f} py={p['y']:.1f}", fill=COL, font=fsm)
    tiles.append(t)
tw, th = tiles[0].size; pad = 10
sheet = Image.new("RGB", (3 * tw + 4 * pad, 3 * th + 4 * pad), (235, 235, 235))
for k, t in enumerate(tiles):
    sheet.paste(t, (pad + (k % 3) * (tw + pad), pad + (k // 3) * (th + pad)))
sheet.save(f"{out}/s35_zoom_tiles.png")
with open(f"{out}/s35_checklist.csv", "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\r\n")
    w.writerow(["no", "px", "py", "k_um_inv", "omega_over_2pi_kHz", "decision", "reason"])
    for i, r in enumerate(rows, 1):
        w.writerow([i, r[2], r[3], r[0], r[1], "", ""])
