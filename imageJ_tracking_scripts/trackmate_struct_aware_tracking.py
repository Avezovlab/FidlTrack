from fiji.plugin.trackmate import TrackMate, Model, Settings, Logger
from fiji.plugin.trackmate.features.edges import EdgeAmbiguityAnalyzer
from fiji.plugin.trackmate.tracking.jaqaman import SimpleSparseLAPTrackerFactory
from fiji.plugin.trackmate.tracking.jaqaman.costfunction import ComponentDistancesTime, ReachableDistCostFunctionTime

from ij import IJ, WindowManager
from fiji.plugin.trackmate import Spot

import os
import sys
from os import path
from math import sqrt


sys.path.append(" path to folder containing a config_tracking.py file")
#eg.
#sys.path.append("../FidlTrack_example_data/240130_cos418+716_3.5ul_6ms")

from config_tracking import *

def load_spots_trackmate(f, model):
	head = f.readline().rstrip("\n").split(",")
	frame_idx = head.index("FRAME")
	x_idx = head.index("POSITION_X")                                 
	y_idx = head.index("POSITION_Y")
	r_idx = head.index("RADIUS")
	q_idx = head.index("QUALITY")
	for i, ln in enumerate(f.readlines()):
		ln = ln.rstrip("\n").split(",")
		frame = int(float(ln[frame_idx]))
		spt = Spot(float(ln[x_idx]), float(ln[y_idx]), 0.0, float(ln[r_idx]), float(ln[q_idx]), "ID{}".format(i))
		spt.putFeature("POSITION_T", float(ln[0]))
		model.addSpotTo(spt, frame)

todo_dirs = []
for exp_dir in os.listdir(out_dir):
	if path.isdir("/".join([out_dir, exp_dir])) and not any([e in exp_dir for e in exclude]):
		todo_dirs.append(exp_dir)

cd = None
cost_f = None

for cpt, exp_path in enumerate(todo_dirs):
	exp_path = "/".join([out_dir, exp_path])
	base_fname = exp_path.split("/")[-1]
	print("Processing[{}/{}]: {}".format(cpt + 1, len(todo_dirs), base_fname))

	comp_f = path.join(base_dir, comp_fname_win.format(fname=base_fname))
	if not path.isfile(comp_f):
		print("   ERROR file not found: {}".format(comp_f))
		continue

	IJ.open(comp_f)
	imp = IJ.getImage()
	dims = imp.getDimensions()
	is_single_frame = all([e == 1 for e in dims[2:]])

	cur_dist_fname = path.join(base_dir, dist_fname.format(fname=base_fname))
	print(cur_dist_fname)
	if not path.isfile(cur_dist_fname):
		print("  Skipped: no distance file found")
		imp.changes = False
		imp.close()
		continue

	settings = Settings()
	settings.tstart = 0
	settings.tend = dims[3]
	settings.dx = set_dx
	settings.dy = set_dx
	settings.dt = set_dt
	settings.width = dims[0]
	settings.height = dims[1]

	settings.addEdgeAnalyzer(EdgeAmbiguityAnalyzer())

	settings.trackerFactory = SimpleSparseLAPTrackerFactory()
	settings.trackerSettings = settings.trackerFactory.getDefaultSettings()


	print("Loading distances")
	cd = ComponentDistancesTime(cur_dist_fname, imp, settings.width, settings.height, settings.dx)
	cost_f = ReachableDistCostFunctionTime(cd)

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

			frames = set()
			cnt = 0
			spots_to_rm = []
			for s in model.getSpots().iterator(True):
				fr = int(s.getFeature("FRAME"))
				frames.add(fr)
				px_pos = [int(s.getFeature("POSITION_X") // set_dx), int(s.getFeature("POSITION_Y") // set_dx)]

				if cd.getComponent(cd.getWindowIdxs(fr)[0], px_pos) == 0:
					spots_to_rm.append(s)
				else:
					s.putFeature("PX_X", px_pos[0])
					s.putFeature("PX_Y", px_pos[1])
				cnt += 1

			print("REMOVING {} over {} spots".format(len(spots_to_rm), cnt))
			for s in spots_to_rm:
				model.removeSpot(s)

			cd.preprocess_spots(model.getSpots())

			#ADD fake spots so that the tracking does not skip inexistant frames
			for i in range(min(frames), max(frames)):
				if i not in frames:
					spt = Spot(float("nan"), float("nan"), 0.0, float("nan"), float("nan"), "ID{}".format(i))
					spt.putFeature("POSITION_T", float("nan"))
					spt.putFeature("PX_X", float("nan"))
					spt.putFeature("PX_Y", float("nan"))
					model.addSpotTo(spt, i)

			settings.tend = model.getSpots().lastKey()
			settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = 0.0
			settings.trackerSettings['MAX_FRAME_GAP'] = 0
			settings.trackerSettings['ALLOW_GAP_CLOSING'] = False

			print(settings.dx, settings.dy, settings.dt, settings.tstart, settings.tend, settings.width, settings.height)

			#======= TRACKING
			for p_DISTANCE in p_LINKING_MAX_DISTANCES:
				tmp = spot_fname.format(mask=mask, spot_rad=p_DIAMETER, spot_th=spot_th)
				if tmp.endswith(".csv"):
					tmp = tmp[:-len(".csv")]
				outFname = path.join(exp_path, link_struct_fname.format(fname=tmp, dist=p_DISTANCE, distgap=0.0, framegap=0))

				if not force and path.isfile(outFname):
					print("  Skipped")
					continue

				print("  Processing: {}".format(outFname))

				# Configure tracker
				settings.trackerSettings['LINKING_MAX_DISTANCE'] = p_DISTANCE
				settings.trackerSettings["COMPONENTS_DISTANCES"] = cd

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
				nerrs = 0
				with open(outFname, 'w') as f:
					f.write('Traj. id, Spot id, x (um), y (um), time (sec), frame, winIdx, component, d_g (um), Ambiguity\n')
					for tid in model.getTrackModel().trackIDs(True):
					    spots = sorted(model.getTrackModel().trackSpots(tid), key=lambda s: s.getFeature('FRAME'))
					    edges = sorted(model.getTrackModel().trackEdges(tid), key=lambda e: model.getTrackModel().getEdgeSource(e).getFeature('FRAME'))

					    for i in range(len(spots)):
					    	spot = spots[i]
					    	if i < len(spots) - 1:
					    		dg = sqrt(cost_f.linkingCost(spots[i], spots[i+1]))
					    		if dg > 1.1 * p_DISTANCE:
					    			nerrs += 1
					    		na = int(model.getFeatureModel().getEdgeFeature(edges[i], "AMBIGUITY"))
				    		else:
				    			dg = -1.0
				    			na = -2

					        sid = spot.ID()
					        x = spot.getFeature('POSITION_X')
					        y = spot.getFeature('POSITION_Y')
					        t = spot.getFeature('POSITION_T')
					        fr =  int(spot.getFeature('FRAME'))

					        s_wins = cd.getWindowIdxs(fr)
					        comp = cd.getComponent(s_wins[0], [int(spot.getFeature("PX_X")), int(spot.getFeature("PX_Y"))])

					        f.write(",".join([str(e) for e in [tid, sid, x, y, t, fr, s_wins[0], comp, dg, na]]) + "\n")
				if nerrs > 0:
					print("  WARNING: {} invalid distances computed, result might be bogus, verify your distance file".format(nerrs))

				model.clearTracks(True)
			model.clearSpots(True)

	imp.changes = False
	imp.close()

	cd = None
	cost_f = None


print("done")
