'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_results.py

讀取 result.csv (row,col,phi_out,l_out,dx,dy,dz)
以及 sky_in.jpg / sky_out.jpg，
根據 l_out>0 決定取樣哪張全景圖，
把方向向量投回 equirectangular 空間，
生成最終渲染圖 final_render.png。
"""

import numpy as np
from PIL import Image
import os

# 1. 參數設定
CSV_PATH     = "result.csv"
SKY_IN_PATH  = "InterstellarWormhole_Fig6a.jpg"
SKY_OUT_PATH = "Cloudysky.jpg"
OUTPUT_PATH  = "final_render.png"

# 2. 載入 CSV
data = np.loadtxt(CSV_PATH, delimiter=',', skiprows=1)
rows   = data[:,0].astype(int)
cols   = data[:,1].astype(int)
l_out  = data[:,3]
dx     = data[:,4]
dy     = data[:,5]
dz     = data[:,6]

# 3. 推斷輸出影像尺寸
height = rows.max() + 1
width  = cols.max() + 1

# 4. 載入全景貼圖
sky1 = Image.open(SKY_IN_PATH).convert('RGB')
sky2 = Image.open(SKY_OUT_PATH).convert('RGB')
arr1 = np.array(sky1, dtype=np.float32) / 255.0  # normalize to [0,1]
arr2 = np.array(sky2, dtype=np.float32) / 255.0
h1, w1, _ = arr1.shape
h2, w2, _ = arr2.shape

# Bilinear sampling function
def bilinear_sample(image, u, v):
    """
    Perform bilinear sampling on image at floating coords u, v.
    image: HxWx3 array in [0,1]
    u: horizontal pixel coordinate (float)
    v: vertical pixel coordinate (float)
    """
    H, W, C = image.shape
    # clamp coordinates
    u = np.clip(u, 0, W - 1)
    v = np.clip(v, 0, H - 1)
    # integer bounds
    x0, y0 = int(np.floor(u)), int(np.floor(v))
    x1, y1 = min(x0 + 1, W - 1), min(y0 + 1, H - 1)
    # fractional parts
    dx, dy = u - x0, v - y0
    # values at corners
    top_left     = image[y0, x0]
    top_right    = image[y0, x1]
    bottom_left  = image[y1, x0]
    bottom_right = image[y1, x1]
    # interpolate
    top    = (1 - dx) * top_left  + dx * top_right
    bottom = (1 - dx) * bottom_left + dx * bottom_right
    return (1 - dy) * top + dy * bottom

# 5. 建立空影像 buffer (float in [0,1])
result = np.zeros((height, width, 3), dtype=np.float32)

# 6. 對每個像素做取樣
for idx in range(data.shape[0]):
    r = rows[idx]
    c = cols[idx]
    vx, vy, vz = dx[idx], dy[idx], dz[idx]
    # 正規化方向向量
    #norm = np.sqrt(vx*vx + vy*vy + vz*vz)
    #if norm == 0:
    #    continue
    #vx /= norm; vy /= norm; vz /= norm

    # equirectangular mapping
    phi = np.arctan2(vy, vx)            # [-π,π]
    u_frac = (phi + np.pi) / (2 * np.pi)  # [0,1)
    v_frac = (vz + 1) * 0.5               # [0,1]

    if l_out[idx] > 0:
        # map to pixel coords in sky1
        u_pix = u_frac * (w1 - 1)
        v_pix = v_frac * (h1 - 1)
        color = bilinear_sample(arr1, u_pix, v_pix)
    else:
        u_pix = u_frac * (w2 - 1)
        v_pix = v_frac * (h2 - 1)
        color = bilinear_sample(arr2, u_pix, v_pix)

    result[r, c] = color

# 7. 存檔 (convert back to uint8)
out_img = (np.clip(result, 0, 1) * 255).astype(np.uint8)
Image.fromarray(out_img).save(OUTPUT_PATH)
print(f"Saved final image to {OUTPUT_PATH}")

'''

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_results.py

讀取 result.csv (row,col,phi_out,l_out,dx,dy,dz)
以及 sky_in.jpg / sky_out.jpg，
根據 l_out>0 決定取樣哪張全景圖，
把方向向量投回 equirectangular 空間，
生成最終渲染圖 final_render.png。
"""

import numpy as np
from PIL import Image
import os

# 1. 參數設定
CSV_PATH     = "result.csv"
SKY_IN_PATH  = "InterstellarWormhole_Fig10.jpg"
SKY_OUT_PATH = "InterstellarWormhole_Fig6a.jpg"
OUTPUT_PATH  = "final_render.png"
# 2. 載入 CSV
#    skiprows=1 跳過表頭
data = np.loadtxt(CSV_PATH, delimiter=',', skiprows=1)
rows     = data[:,0].astype(int)
cols     = data[:,1].astype(int)
l_out    = data[:,2]
theta_out= data[:,3]
phi_out  = data[:,4]

# 3. 推斷輸出影像尺寸
height = rows.max() + 1
width  = cols.max() + 1

# 4. 載入全景貼圖
sky1 = Image.open(SKY_IN_PATH).convert('RGB')
sky2 = Image.open(SKY_OUT_PATH).convert('RGB')
arr1 = np.array(sky1)
arr2 = np.array(sky2)
h1, w1, _ = arr1.shape
h2, w2, _ = arr2.shape

# 5. 建立空影像 buffer
result = np.zeros((height, width, 3), dtype=np.uint8)

# 6. 對每個像素做取樣（nearest‐neighbor）
for idx in range(data.shape[0]):
    r     = rows[idx]
    c     = cols[idx]
    theta = theta_out[idx]   # ∈ [0, π]
    phi   = phi_out[idx]     # ∈ [−π, π]

    # equirectangular → UV
    u_frac = (phi + np.pi) / (2 * np.pi)  # -> [0, 1)
    v_frac = theta / np.pi                # -> [0, 1]

    if l_out[idx] > 0:
        # 前方全景圖
        iu = int(u_frac * (w1 - 1)) % w1
        iv = int(v_frac * (h1 - 1)) % h1
        color = arr1[iv, iu]
    else:
        # 背方全景圖
        iu = int(u_frac * (w2 - 1)) % w2
        iv = int(v_frac * (h2 - 1)) % h2
        color = arr2[iv, iu]

    result[r, c] = color
# 7. 存檔
out_img = Image.fromarray(result)
out_img.save(OUTPUT_PATH)
print(f"Saved final image to {OUTPUT_PATH}")
