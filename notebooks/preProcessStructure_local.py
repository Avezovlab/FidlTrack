#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  7 12:45:05 2025

@author: pierre
"""
from os import path
import numpy as np
import skimage as sk
from math import floor
from skimage.measure import label, regionprops

img_dir = "XX" #path to the directory containing the data
img_name = "XX" #name of the segmented image


pxsize = XX #in micron
stab_Nframes = XX #in frames
w_dur = XX #in frames
w_ovlp = XX #in frames
rm_comps_ltpxs = XX #in pixels


force_recompute = False


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