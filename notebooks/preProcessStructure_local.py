#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 12:45:05 2025

@author: pierre
"""
from os import path
import numpy as np
import skimage as sk
import struct
from math import floor
import progressbar
from skimage.measure import label, regionprops

import queue
from dataclasses import dataclass, field
from typing import Any

from matplotlib import pyplot as plt

#img_dir = "/mnt/data4/SPT_method_moved_for_space/roger/Hela_250226/CTL"
#img_name = "C2-250226_HeLa_Sec13_SNAP_GFP_Sec61_Halo_KDEL_50nM_PAJF646_c3.nd2_preview_Simple_Segmentation_bin.tif"

#img_dir = "/mnt/data4/SPT_method_moved_for_space/roger/Hela_250226/CTL"
#img_name = "C2-250226_HeLa_Sec13_SNAP_GFP_Sec61_Halo_KDEL_250nM_PAJF646_c9.nd2_preview_Simple Segmentation_bin_open1pxcirc.tif"

# max_dist = 1.3
# w_dur = 21
# w_ovlp = 0.04
# stab_Nframes = 1
# rm_comps_ltpxs = 50
# pxsize = 0.0967821
# force_recompute = False


#img_dir = "/mnt/data2/SPT_method/simu/hex/sim"
#img_name = "hexnet_25_100_poly.poly_fov_dil.tif" #"hexnet_75_100_poly.poly_fov_dil.tif" #"hexnet_38_100_poly.poly_fov_dil.tif"

#img_dir = "/mnt/data2/SPT_method/simu/hex_deci/"
#img_name = "hexnet_deci_25_100_0.025_0.15_bin_fov_dil.tif" #"hexnet_75_100_poly.poly_fov_dil.tif" #"hexnet_38_100_poly.poly_fov_dil.tif"


# img_dir = "/mnt/data2/SPT_method/simu/hex_deci/a"
# img_name = "hexnet_deci_100_100_0.025_0.15_bin_fov_dil.tif"
# pxsize = 0.006048881


#img_dir = "/mnt/data4/SPT_method_moved_for_space/yutong_240123/240123_Yutong_dATL_20ms/cell6/sim"
#img_name = "C1-cell6_MMStack_Pos0_c.ome.tif_avg51_FRAME2252_usharp2px_0.8_blur0.5px_Simple_Segmentation_bin_erodecric1px_adj_dil.tif"
#pxsize = 0.0645

#img_dir = "/mnt/data4/SPT_method_moved_for_space/yutong_240123/240122_Yutong_cos123-716-717-418_HPA646-3ul_6ms/cell1"
#img_name = "C1-cell1_MMStack_Pos0_c.ome.tif_avg51_Simple_Segmentation_binary_cleaned2_closed_eroded_circ1px_inv.tif"

#img_dir = "/mnt/data4/SPT_method_moved_for_space/APP/290725_PP_YY_cos123-931/C1-cell5_2_10ms"
#img_name = "C1-cell5_2_10ms_MMStack_Pos0.ome.tif_avg51_Simple_Segmentation_bin_closed_eroded_circ1px.tif"

#img_dir = "/mnt/data4/SPT_method_moved_for_space/APP/290725_PP_YY_cos123-931/cell10_nobace1_10ms"
#img_name = "C1-cell10_nobace1_10ms_MMStack_Pos0.ome.tif_avg51_Simple_Segmentation_close_erode_circ1px_cleaned.tif"
#pxsize = 0.0645

img_dir = "/mnt/data4/SPT_method_moved_for_space/APP/290725_PP_YY_cos123-931/cell11_nobace1_1_10ms"
img_name = "C1-cell11_nobace1_1_10ms_MMStack_Pos0.ome.tif_avg51_Simple_Segmentation_bin_closed_circ1px.tif" #rm_comps_ltpxs = 200
pxsize = 0.065

stab_Nframes = 3
w_dur = 101
w_ovlp = 0.0001
#w_dur = 21
#w_ovlp = 0.0005
rm_comps_ltpxs = 200
force_recompute = False


#img_dir = "/mnt/data2/SPT_method/simu/lines/"
#img_name = "struct_line_dist=52_pxsize=0.024195525_poly_fov.tif" #"struct_line_dist=42_pxsize=0.024195525_poly_fov.tif" #"struct_line_dist=31_pxsize=0.024195525_poly_fov.tif"


# w_dur = 60001
# w_ovlp = 0
# stab_Nframes = 1
# rm_comps_ltpxs = 0
# #pxsize = 0.024195525
# force_recompute = False


mask_img = "/".join([img_dir, img_name])


if w_ovlp >= 1:
  print("ERROR 0 <= w_ovlp < 1")
  assert(False)

dw = floor(w_dur * (1 - w_ovlp))
if w_ovlp == 0:
    dw = w_dur
print("Window length = {} frame(s)".format(dw))
print("Overlap = {} frame(s)".format(w_dur - dw))


stab_fname = path.join(img_dir, "{}_stabN={}".format(path.splitext(img_name)[0], stab_Nframes))
comps_fname = "{}_comps".format(stab_fname)
win_fname = "{}_wDur={}_wOvlp={}".format(comps_fname, w_dur, w_ovlp)

if not path.isfile(mask_img):
  print("ERROR mask_img file not found: {}".format(mask_img))
  assert(False)



if not force_recompute and path.isfile(stab_fname + ".tif"):
  img = sk.io.imread(stab_fname + ".tif")
else:
  img = sk.io.imread(mask_img)
  if img.ndim == 2:
    img = img.reshape((1, img.shape[0], img.shape[1]))

  if stab_Nframes > 1:
    res = np.zeros(img.shape, img.dtype)
    for i in range(img.shape[0]):
      n_left = int(stab_Nframes/2)
      if i < stab_Nframes/2:
          n_left = i
      n_right = int(stab_Nframes/2) + 1
      if i >= img.shape[0] - stab_Nframes/2:
          n_right = img.shape[0] - i

      cur_imgs = img[(i-n_left):(i+n_right), :, :]
      res[i,:,:] = np.max(cur_imgs, axis=0)
    img = res

  sk.io.imsave(stab_fname + ".tif", img, check_contrast=False)
  print("Saved: {}".format(stab_fname + ".tif"))

print("Raw stack length = {}".format(img.shape[0]))


if not force_recompute and path.isfile(win_fname + ".tif"):
  img = sk.io.imread(win_fname + ".tif")
else:
  win_start = [0] + list(range(dw, img.shape[0], dw))

  labs = np.zeros((len(win_start), img.shape[1], img.shape[2]),
                  dtype="uint16")
  for k, wc in enumerate(win_start):
    for i in range(wc, min(wc + w_dur, img.shape[0])):
      labs[k] = np.maximum(labs[k], img[i] > 0)

    labs[k] = label(labs[k] > 0)
    idxs = np.array(sorted(set(labs[k].flatten()), key=lambda x: sum(labs[k].flatten() == x),
                  reverse=True))

    for p in regionprops(labs[k]):
      if p.coords.shape[0] < rm_comps_ltpxs:
        labs[k, p.coords[:,0], p.coords[:,1]] = 0
      else:
        print(p.label, np.where(idxs == p.label)[0][0])
        labs[k, p.coords[:,0], p.coords[:,1]] = np.where(idxs == p.label)[0][0] + 1

  sk.io.imsave(win_fname + ".tif", labs, check_contrast=False)
  print("Saved: {}".format(win_fname + ".tif"))
  with open(win_fname + "_winstart.csv", 'w') as f:
    f.write(",".join([str(e) for e in win_start]) + "\n")
  print("Saved: {}".format(win_fname + "_winstart.csv"))

print("2D mask? {}".format("yes" if labs.shape[0] == 1 else "no"))
print("Number of time windows: {}".format(labs.shape[0]))