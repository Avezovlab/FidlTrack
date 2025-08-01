from fiji.plugin.trackmate import TrackMate, Model, Settings, Logger
from fiji.plugin.trackmate.tracking.jaqaman import SimpleSparseLAPTrackerFactory
from fiji.plugin.trackmate.features.edges import EdgeAmbiguityAnalyzer

from ij import IJ
from fiji.plugin.trackmate import Spot

import os
import sys
from os import path

#sys.path.append("/mnt/data4/yutong/2cols/240816_cos123-663-681_2colorSPT")
#sys.path.append("/mnt/data4/yutong/2cols/240816_cos123-663-681_1reglocPA")
#sys.path.append("/mnt/data4/yutong/2cols/240816_cos123-663-600_2colorSPT")
#sys.path.append("/mnt/data4/yutong/2cols/240816_cos123-663-600_1reglocPA")
#sys.path.append("/mnt/data3/droso/data/271023_droso/from_241023/SPT")

#sys.path.append("/mnt/data3/droso/data/271023_droso/from_241023/SPT")
#sys.path.append("/mnt/data3/droso/data/140224_PierreMaysoon_drosoMutant_d+2/good")
#sys.path.append("/mnt/data3/droso/data/310124_PierreMaysoon_mutant/good")
#sys.path.append("/mnt/data3/droso/data/010324_droso_3xmutant/good")
#sys.path.append("/mnt/data3/droso/data/201023_droso_cool/DT0.012")
#sys.path.append("/mnt/data3/droso/data/201023_droso_cool/DT0.0135")
#sys.path.append("/mnt/data3/droso/data/pierre_maysoon_280624_droso/wt_3/SPT")
#sys.path.append("/mnt/data3/droso/data/pierre_maysoon_droso_WT_030724/wt1/SPT")

#sys.path.append("/mnt/data3/perk_ire1/raw_data/DRI/010523_U2OS25")
#sys.path.append("/mnt/data3/perk_ire1/raw_data/DRI/040423_H3")

#sys.path.append("/mnt/data4/yutong/2cols/240823_cos123-663-681_2colorSPT/6ms")
#sys.path.append("/mnt/data4/yutong/2cols/240823_cos123-663-681_2colorSPT/12ms")
#sys.path.append("/mnt/data4/yutong/2cols/240823_cos123-663-600_2colorSPT/6ms")


#sys.path.append("/mnt/data4/SPT_ineurons/data/230717_Yutong_neuron10_SPT")
#sys.path.append("/mnt/data4/SPT_ineurons/data/230714_Yutong_neuron7_SPT")
#sys.path.append("/mnt/data4/SPT_ineurons/data/230808_Yutong_neuron13_SPT_6ms")
#sys.path.append("/mnt/data4/SPT_ineurons/data/230818_Yutong_neuron_SPT_multiregion")

#sys.path.append("/home/pierre/yutong/bip/051124")
#sys.path.append("/home/pierre/yutong/bip/191124_BiP_2colorSPT")

#sys.path.append("/mnt/data2/SPT_method/nanobody/869_717/250116")
#sys.path.append("/mnt/data2/SPT_method/roger/u2os_HaloKDEL_250123")
#sys.path.append("/mnt/data2/SPT_method/nanobody/869_717/250120")


#sys.path.append("/mnt/data2/SPT_method/roger/Hela_250206")
#sys.path.append("/mnt/data4/SPT_method_moved_for_space/roger/Hela_250226")
#sys.path.append("//mnt/data2/SPT_method/roger/Hela_250220")

#sys.path.append("/mnt/data2/SPT_method/roger/u2os_HaloKDEL")

#sys.path.append("/mnt/data2/SPT_method/nanobody/APP")
#sys.path.append("/mnt/data2/SPT_method/nanobody/nb")
sys.path.append("/mnt/data2/SPT_method/nanobody/nb+APP")

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
for root, dirs, files in os.walk(out_dir):
	if not any([e in root for e in exclude]) and base_start in root and any([e.startswith("spots_") for e in files]): 
		todo_dirs.append(root)

print(base_dir)
for cpt, exp_path in enumerate(todo_dirs):
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
		if type(ths) == dict:
			cur_ths = ths[base_fname[:2]]
		else:
			cur_ths = ths

		for spot_th in cur_ths:
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
					
					print("  Processing: {} dist={} gap={}".format(tmp, p_DISTANCE, p_GAP_FRAME))

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
