#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 14:17:32 2025

@author: pierre
"""
import numpy as np
import struct
import progressbar
from skimage.measure import regionprops
from skimage.io import imread

import queue
from dataclasses import dataclass, field
from typing import Any


@dataclass(order=True)
class Item:
    priority: int
    v: Any=field(compare=False)

def d1D(px1D, ppx1D, IMSIZE, SQRT2_PXSIZE, PXSIZE):
    return SQRT2_PXSIZE if abs(px1D - ppx1D) in {IMSIZE-1, IMSIZE+1} else PXSIZE

def neighbors_1D(px1D, IMSIZE):
  px = [px1D // IMSIZE, px1D % IMSIZE]
  nh = [[px[0] + nh[0], px[1] + nh[1]] for nh in [[0, -1], [-1, 0], [+1, 0], [0, +1],
                                                  [-1, -1], [+1, -1], [-1, +1], [+1, +1]]]
  return [e[0] * IMSIZE + e[1] for e in nh if 0 <= e[0] <= (IMSIZE-1) and 0 <= e[1] <= (IMSIZE-1)]

def time_ranges_ovlp(tr1, tr2):
  min_ovlp = max(tr1[0], tr2[0])
  max_ovlp = min(tr1[1], tr2[1])
  if max_ovlp >= min_ovlp:
    return [min_ovlp, max_ovlp]
  return None

def find_neighbors3D(px1D, comp, ts, te, pxs_grps, pxs3D, pxs1D_to_3D, imsize):
  res = []
  for nh_px1D in neighbors_1D(px1D, imsize):
    if nh_px1D not in pxs1D_to_3D:
      continue
    for p3D in pxs1D_to_3D[nh_px1D]:
      p1D,j = pxs3D[p3D]
      assert(p1D == nh_px1D)
      nh_ts,nh_te,nh_c = pxs_grps[p1D][j]
      t_ovlp = time_ranges_ovlp([ts, te], [nh_ts, nh_te])
      if nh_c == comp and t_ovlp != None:
        res.append((p3D, p1D, comp, t_ovlp[0], t_ovlp[1]))
  return res

def iscontained(tr1, tr2):
  return (tr1[0] >= tr2[0] and tr1[0] <= tr2[1] and
          tr1[1] <= tr2[1] and tr1[1] >= tr2[0])

def split_timerange(tr1, tr2):
  res = []
  min_ovlp = max(tr1[0], tr2[0])
  max_ovlp = min(tr1[1], tr2[1])

  if max_ovlp < min_ovlp:
    return []

  if tr1[0] < min_ovlp:
    res.append([1, tr1[0], min_ovlp-1])
  elif tr2[0] < min_ovlp:
    res.append([2, tr2[0], min_ovlp-1])
  res.append([12, min_ovlp, max_ovlp])
  if tr1[1] > max_ovlp:
    res.append([1, max_ovlp+1, tr1[1]])
  elif tr2[1] > max_ovlp:
    res.append([2, max_ovlp+1, tr2[1]])
  return res

def merge_ranges(new_elts):
  if len(new_elts) == 1:
    return new_elts

  merge_done = False
  while not merge_done:
    merge_done = True
    to_del = []
    for u in range(len(new_elts)):
      tr1 = new_elts[u][1:3]
      d1 = new_elts[u][0]
      for v in range(len(new_elts)):
        if u == v:
          continue
        tr2 = new_elts[v][1:3]
        d2 = new_elts[v][0]

        if iscontained(tr1, tr2) and d2 <= d1:
          to_del = [u]
          merge_done = False
          continue
        if iscontained(tr2, tr1) and d1 <= d2:
          to_del = [v]
          merge_done = False
          continue

        isets = split_timerange(tr1, tr2)
        if isets == []:
          continue
        to_del = [max([u,v]), min([u,v])]
        merge_done = False
        for iset in isets:
          if iset[0] == 1:
            new_elts.append([d1, *iset[1:3]])
          elif iset[0] == 2:
            new_elts.append([d2, *iset[1:3]])
          elif iset[0] == 12:
            new_elts.append([min([d1, d2]), *iset[1:3]])
          else:
            assert(False)
      if not merge_done:
        break
    for v in to_del:
      del new_elts[v]

  return new_elts

def dijkstra_2D(px1D, pxs1D, max_dist, imsize, pxsize):
  s2pxsize = np.sqrt(2) * pxsize
  idxs = {pxs1D[i]: i for i in range(len(pxs1D))}
  nh1D = np.array([-imsize, -1, +1, imsize, imsize-1, imsize+1,
                   -imsize-1, -imsize+1])

  dists = np.ones((len(pxs1D))) * (2 * max_dist)
  dists[idxs[px1D]] = 0
  q = queue.PriorityQueue()
  q.put(Item(dists[idxs[px1D]], px1D))
  done = np.zeros(len(pxs1D), dtype="bool")

  spxs1D = set(pxs1D)

  while not q.empty():
    cur_px1D = q.get().v

    if done[idxs[cur_px1D]]:
      continue

    done[idxs[cur_px1D]] = 1
    cur_px2D = [cur_px1D // imsize, cur_px1D % imsize]

    for nh_px1D in [e for e in cur_px1D+nh1D if e in spxs1D and done[idxs[e]] == 0]:
      #check if we loop through the image
      nh_px2D = [nh_px1D // imsize, nh_px1D % imsize]
      if abs(cur_px2D[0] - nh_px2D[0]) > 1 or abs(cur_px2D[1] - nh_px2D[1]) > 1:
        continue

      nh_d = dists[idxs[cur_px1D]] + d1D(cur_px1D, nh_px1D, imsize, s2pxsize, pxsize)

      if nh_d < max_dist and nh_d < dists[idxs[nh_px1D]]:
        dists[idxs[nh_px1D]] = nh_d
        q.put(Item(nh_d, nh_px1D))
  return dists

def save_bin_comps2D(all_dists, cache_r, dists, rev_map_dists, out_fname):
  with open(out_fname, 'wb') as f:
    f.write(int(2).to_bytes(3, byteorder='big'))

    f.write(len(all_dists).to_bytes(3, byteorder='big'))
    for d in all_dists:
      f.write(struct.pack('>f', d))

    f.write(len(dists).to_bytes(3, byteorder='big'))
    for lab in sorted(dists.keys()):
      f.write(lab.to_bytes(3, byteorder='big'))
      f.write(len(dists[lab]).to_bytes(3, byteorder='big'))
      for px3D1 in sorted(dists[lab].keys()):
        f.write(int(px3D1).to_bytes(3, byteorder='big'))
        f.write(len(dists[lab][px3D1]).to_bytes(3, byteorder='big'))
        for px3D2 in sorted(dists[lab][px3D1].keys()):
          f.write(int(px3D2).to_bytes(3, byteorder='big'))
          f.write(rev_map_dists[cache_r[dists[lab][px3D1][px3D2]]].to_bytes(3, byteorder='big'))

def save_bin_comps3D(w_dur, w_ovlp, nframes, all_dists, cache_r, dists, rev_map_dists, out_fname):
  with open(out_fname, 'wb') as f:
    f.write(int(1).to_bytes(3, byteorder='big'))
    f.write(w_dur.to_bytes(3, byteorder='big'))
    f.write(struct.pack('>f', w_ovlp))
    f.write(int(nframes).to_bytes(3, byteorder='big'))

    f.write(len(all_dists).to_bytes(3, byteorder='big'))
    for d in all_dists:
      f.write(struct.pack('>f', d))

    f.write(len(dists).to_bytes(3, byteorder='big'))
    for lab in sorted(dists.keys()):
      f.write(lab.to_bytes(3, byteorder='big'))
      f.write(len(dists[lab]).to_bytes(3, byteorder='big'))
      for px1 in sorted(dists[lab].keys()):
        f.write(int(px1).to_bytes(3, byteorder='big'))
        f.write(len(dists[lab][px1]).to_bytes(3, byteorder='big'))
        for px2 in sorted(dists[lab][px1].keys()):
          f.write(int(px2).to_bytes(3, byteorder='big'))
          f.write(len(dists[lab][px1][px2]).to_bytes(3, byteorder='big'))
          for elts in dists[lab][px1][px2]:
            f.write(rev_map_dists[cache_r[elts[0]]].to_bytes(3, byteorder='big'))
            f.write(int(elts[1]).to_bytes(3, byteorder='big'))
            f.write(int(elts[2]).to_bytes(3, byteorder='big'))

####STARTS HERE


basedir = "XX" #path to directory containing the data
fname = "XX" #name of the component stack file


pxsize = XX #in micron
max_dist = XX #in micron
stab_Nframes = XX
w_dur = XX #in frames
w_ovlp = XX
rm_comps_ltpxs = XX #in pixels
force_recompute = False

dist_fname = "{}_dist={}".format(fname, max_dist)
SQRT2_PXSIZE = np.sqrt(2) * pxsize


labs = imread("{}/{}".format(basedir, fname))
if len(labs.shape) == 2:
    labs = labs.reshape((1, labs.shape[0], labs.shape[1]))

IMSIZE = labs.shape[2]

print(fname)

elts_f = {}
elts_r = {}
dists = {}
if labs.shape[0] == 1: #2D mask
  all_px1Ds = set()
  comps_px1D = {}
  pr = regionprops(labs[0])
  for p in pr:
    px1Ds = set(p.coords[:, 0] * IMSIZE + p.coords[:, 1])
    all_px1Ds.update(px1Ds)
    comps_px1D[p.label] = set(px1Ds)
  all_px1Ds = sorted(all_px1Ds)

  for k in comps_px1D.keys():
    comps_px1D[k] = sorted(comps_px1D[k])

  for comp,pxs in comps_px1D.items():
    comp = int(comp)
    dists[comp] = {}
    for k in progressbar.progressbar(range(len(pxs))):
      tmp = dijkstra_2D(pxs[k], pxs, max_dist, IMSIZE, pxsize)
      tmp = {pxs[i]: np.round(tmp[i],6) for i in range(len(tmp)) if tmp[i] < max_dist and pxs[i] > pxs[k]}

      for d in tmp.values():
        if d not in elts_f:
          ncache = len(elts_f)+1
          elts_f[d] = ncache
          elts_r[ncache] = d
      dists[comp][pxs[k]] = {p: elts_f[d] for p,d in tmp.items()}

  all_dists = set()
  for comp, elts1 in dists.items():
    for px3D1, elts2 in elts1.items():
      all_dists.update([elts_r[didx] for didx in elts2.values()])

  all_dists = sorted(all_dists)
  rev_map_dists = {}
  for i,d in enumerate(all_dists):
    rev_map_dists[d] = i

  save_bin_comps2D(all_dists, elts_r, dists, rev_map_dists, dist_fname + ".bin")
  print("Saved: {}".format(dist_fname + ".bin"))
else: #stack of masks
  for n in progressbar.progressbar(range(labs.shape[0])):
    all_px1Ds = set()
    comps_px1D = {}
    pr = regionprops(labs[n])
    for p in pr:
      px1Ds = set(p.coords[:, 0] * IMSIZE + p.coords[:, 1])
      all_px1Ds.update(px1Ds)
      comps_px1D[p.label] = set(px1Ds)
    all_px1Ds = sorted(all_px1Ds)

    for k in comps_px1D.keys():
      comps_px1D[k] = sorted(comps_px1D[k])

    for comp,pxs in comps_px1D.items():
      print(" ", len(pxs))
      comp = int(comp)
      if comp not in dists:
        dists[comp] = {}

      for k in range(len(pxs)):
        tmp = dijkstra_2D(pxs[k], pxs, max_dist, IMSIZE, pxsize)
        tmp = {pxs[i]: np.round(tmp[i],6) for i in range(len(tmp)) if tmp[i] < max_dist and pxs[i] > pxs[k]}

        if tmp != {} and pxs[k] not in dists[comp]:
          dists[comp][pxs[k]] = {}

        for p2,d in tmp.items():
          if d not in elts_f:
            ncache = len(elts_f)+1
            elts_f[d] = ncache
            elts_r[ncache] = d

          if p2 not in dists[comp][pxs[k]]:
            dists[comp][pxs[k]][p2] = {}
          if elts_f[d] not in dists[comp][pxs[k]][p2]:
            dists[comp][pxs[k]][p2][elts_f[d]] = []
          dists[comp][pxs[k]][p2][elts_f[d]].append(n)

      if len(dists[comp]) == 0:
        del dists[comp]

  dists_tmp = dists

  for c in dists.keys():
    for p1 in dists[c].keys():
      for p2 in dists[c][p1].keys():
        tmp = []
        for d,elts in dists[c][p1][p2].items():
          start = elts[0]
          elts = sorted(elts)
          for i in range(1, len(elts)):
            if elts[i-1] != elts[i] - 1:
              tmp.append((d, start, elts[i-1]))
              start = elts[i]
          tmp.append((d, start, elts[len(elts)-1]))
        dists[c][p1][p2] = tmp

  all_dists = set()
  for comp, elts1 in dists.items():
    for px3D1, elts2 in elts1.items():
      for px3D2, elts3 in elts2.items():
        [all_dists.add(elts_r[elt[0]]) for elt in elts3]

  all_dists = sorted(all_dists)
  rev_map_dists = {}
  for i,d in enumerate(all_dists):
    rev_map_dists[d] = i

  save_bin_comps3D(w_dur, w_ovlp, labs.shape[0], all_dists, elts_r, dists, rev_map_dists, dist_fname + ".bin")
  print("Saved: {}".format(dist_fname + ".bin"))