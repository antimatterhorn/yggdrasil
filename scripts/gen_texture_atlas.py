import colorsys
import os
from PIL import Image

WIDTH = 1024
HEIGHT = 1024
SWATCH_W = 32
SWATCH_H = 64

num_cols = WIDTH // SWATCH_W        # 32 columns
num_rows = (HEIGHT // 4) // SWATCH_H  # 8 rows (top half only)
total_swatches = num_cols * num_rows   # 256 swatches

img = Image.new("RGBA", (WIDTH, HEIGHT), (0, 0, 0, 0))
pixels = img.load()

HUE_SPEED = num_rows  # each row independently completes one full hue cycle

for idx in range(total_swatches):
    row = idx // num_cols
    col = idx % num_cols
    row_t = row / (num_rows - 1)  # 0.0 at row 0, 1.0 at bottom row

    hue = (idx / total_swatches * HUE_SPEED) % 1.0

    # Progressively desaturate toward the bottom rows (vivid → muted/earthy)
    row_sat_mult = 1.0 - row_t * 0.85   # 1.00 at row 0 → 0.15 at row 7

    # Compress and lower the lightness range toward the bottom rows (darker overall)
    light_max = 0.90 - row_t * 0.40     # 0.90 at row 0 → 0.50 at row 7
    light_min = 0.15 - row_t * 0.07     # 0.15 at row 0 → 0.08 at row 7

    for py in range(SWATCH_H):
        t = 1.0 - py / (SWATCH_H - 1)  # 1=top (light), 0=bottom (dark)
        lightness = light_min + t * (light_max - light_min)
        saturation = (1.0 - t * 0.25) * row_sat_mult
        r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
        color = (int(r * 255), int(g * 255), int(b * 255), 255)
        border = py == 0 or py == SWATCH_H - 1
        for px in range(SWATCH_W):
            is_border = border or px == 0 or px == SWATCH_W - 1
            pixels[col * SWATCH_W + px, row * SWATCH_H + py] = (0, 0, 0, 255) if is_border else color

# Row 5 (index 4): gray reference swatches — saturation = 0
# (light_top, light_bottom): uniform chips have equal values; last is the full gradient
GRAY_ROW = 4
gray_swatches = [
    (1.00, 0.70),  # white → light gray
    (1.00, 0.50),  # white → light gray
    (1.00, 0.30),  # white → light gray
    (1.00, 0.00),  # white → black (full range)
    (0.80, 0.00),  # light gray → medium gray
    (0.60, 0.00),  # medium gray → medium-dark gray
    (0.40, 0.00),  # medium-dark gray → black
]

for col_idx, (light_top, light_bot) in enumerate(gray_swatches):
    for py in range(SWATCH_H):
        t = 1.0 - py / (SWATCH_H - 1)  # 1=top (light), 0=bottom (dark)
        lightness = light_bot + t * (light_top - light_bot)
        lv = int(lightness * 255)
        border = py == 0 or py == SWATCH_H - 1
        for px in range(SWATCH_W):
            is_border = border or px == 0 or px == SWATCH_W - 1
            pixels[col_idx * SWATCH_W + px, GRAY_ROW * SWATCH_H + py] = (0, 0, 0, 255) if is_border else (lv, lv, lv, 255)

out = os.path.join(os.path.dirname(__file__), "texture_atlas.png")
img.save(out)
print(f"Saved {out}  ({num_cols}×{num_rows} swatches, {SWATCH_W}×{SWATCH_H}px each)")
