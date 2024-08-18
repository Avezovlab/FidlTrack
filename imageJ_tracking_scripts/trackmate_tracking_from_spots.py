from fiji.plugin.trackmate import TrackMate, Model, Settings, Logger
from fiji.plugin.trackmate.tracking.jaqaman import SimpleSparseLAPTrackerFactory
from fiji.plugin.trackmate.features.edges import EdgeAmbiguityAnalyzer

from ij import IJ
from fiji.plugin.trackmate import Spot

import os
import sys
from os import path

sys.path.append(" path to folder containing a config_tracking.py file")
#eg.
#sys.path.append("../FidlTrack_example_data/240130_cos418+716_3.5ul_6ms")

from config_tracking import *

def load_spots_trackmate(f, model):
	for i, ln in enumerate(f.readlines()):
		ln = ln.rstrip("\n").split(",")
		if i == 0:
			continue
		frame = int(float(ln[3]))
		spt = Spot(float(ln[1]), float(ln[2]), 0.0, float(ln[4]), float(ln[5]), "ID{}".format(i))
		spt.putFeature("POSITION_T", float(ln[0]))
		model.addSpotTo(spt, frame)

todo_dirs = []
todo_dirs = []
for exp_dir in os.listdir(out_dir):
	if path.isdir("/".join([out_dir, exp_dir])) and not any([e in exp_dir for e in exclude]):
		todo_dirs.append(exp_dir)

print(base_dir)
for cpt, exp_path in enumerate(todo_dirs):
	exp_path = "/".join([out_dir, exp_path])
	base_fname = exp_path.split("/")[-1]
	print("Processing[{}/{}]: {}".format(cpt + 1, len(todo_dirs), base_fname))

	settings = Settings()
	settings.tstart = 0
	settings.dx = set_dx
	settings.dy = set_dx
	settings.dt = set_dt

	settings.addEdgeAnalyzer(EdgeAmbiguityAnalyzer())

	settings.trackerFactory = SimpleSparseLAPTrackerFactory()
	settings.trackerSettings = settings.trackerFactory.getDefaultSettings()

	print("pxsize = {} um; dt = {} s".format(settings.dx, settings.dt))

	for p_DIAMETER in p_DETECTION_BLOB_DIAMETERS:
		for spot_th in ths:
			fname = "/".join([exp_path, spot_fname.format(mask=mask, spot_rad=p_DIAMETER, spot_th=spot_th)])
			if not path.isfile(fname):
				print("  Skipped: spot file not found: {}".format(fname))
				continue
			else:
				print("  {}".format(fname))

			model = Model()
			model.setLogger(Logger.IJ_LOGGER)
			model.setPhysicalUnits("µm", "ms")

			with open(fname, 'r') as f:
				load_spots_trackmate(f, model)

			#ADD fake spots so that the tracking does not skip inexistant frames
			frames = set([int(s.getFeature("FRAME")) for s in model.getSpots().iterator(True)])
			for i in range(min(frames), max(frames)):
				if i not in frames:
					spt = Spot(float("nan"), float("nan"), 0.0, float("nan"), float("nan"), "ID{}".format(i))
					spt.putFeature("POSITION_T", float("nan"))
					spt.putFeature("PX_X", float("nan"))
					spt.putFeature("PX_Y", float("nan"))
					model.addSpotTo(spt, i)

			settings.tend = model.getSpots().lastKey()
			trackmate = TrackMate(model, settings)
			trackmate.setNumThreads(4)
			
			print(settings.dx, settings.dy, settings.dt, settings.tstart, settings.tend, settings.width, settings.height)

			#======= TRACKING
			for p_DISTANCE in p_LINKING_MAX_DISTANCES:
				for p_GAP_FRAME in p_MAX_FRAME_GAPS:
					if p_GAP_FRAME == 0:
						gap_distance = 0.0
					else:
						gap_distance = p_DISTANCE

					tmp = spot_fname.format(mask=mask, spot_rad=p_DIAMETER, spot_th=spot_th)
					if tmp.endswith(".csv"):
						tmp = tmp[:-len(".csv")]
					outFname = "/".join([exp_path, link_fname.format(fname=tmp, dist=p_DISTANCE, distgap=gap_distance, framegap=p_GAP_FRAME)])

					if not force and path.isfile(outFname):
						print("  Skipped")
						continue
					
					print("  Processing: {}".format(tmp))

					# Configure tracker
					settings.trackerSettings['LINKING_MAX_DISTANCE'] = p_DISTANCE
					settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = gap_distance
					settings.trackerSettings['MAX_FRAME_GAP'] = p_GAP_FRAME

					if p_GAP_FRAME == 0:
						settings.trackerSettings['ALLOW_GAP_CLOSING'] = False

					trackmate = TrackMate(model, settings)
					trackmate.setNumThreads(4)

					ok = trackmate.execTracking()
					if not ok:
						print(str(trackmate.getErrorMessage()))
						continue
					ok = trackmate.computeTrackFeatures(True)
					if not ok:
						print(str(trackmate.getErrorMessage()))
						continue
					ok = trackmate.execTrackFiltering(True)
					if not ok:
						print(str(trackmate.getErrorMessage()))
						continue
					ok = trackmate.computeEdgeFeatures(True)
					if not ok:
						print(str(trackmate.getErrorMessage()))
						continue

					#================= EXPORT TRACKS
					fm = model.getFeatureModel()
					with open(outFname, 'w') as f:
						f.write('Traj. id, Spot id, x (um), y (um), time (sec), frame, Ambiguity\n')
						for tid in model.getTrackModel().trackIDs(True):
						    spots = sorted(model.getTrackModel().trackSpots(tid), key=lambda s: s.getFeature('FRAME'))
						    edges = sorted(model.getTrackModel().trackEdges(tid), key=lambda e: model.getTrackModel().getEdgeSource(e).getFeature('FRAME'))

						    for i in range(len(spots)):
						    	spot = spots[i]

						    	if i < len(spots) - 1:
						    		na = int(model.getFeatureModel().getEdgeFeature( edges[i], "AMBIGUITY" ))
					    		else:
					    			na = -2

						        sid = spot.ID()
						        x = spot.getFeature('POSITION_X')
						        y = spot.getFeature('POSITION_Y')
						        t = spot.getFeature('POSITION_T')
						        fr = spot.getFeature('FRAME')

						        f.write(",".join([str(e) for e in [tid, sid, x, y, t, fr, na]]) + "\n")

					model.clearTracks(True)
			model.clearSpots(True)

print("done")
