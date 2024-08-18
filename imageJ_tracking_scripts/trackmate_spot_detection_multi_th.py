from fiji.plugin.trackmate import TrackMate, Settings, Model, Logger, Spot
from fiji.plugin.trackmate.detection import LogDetectorFactory

from ij import IJ, WindowManager
import fiji.plugin.trackmate.features.FeatureFilter as FeatureFilter

from ij.measure import Calibration

import os
import sys
from os import path
from math import sqrt

sys.path.append(" path to folder containing a config_tracking.py file")
#eg.
#sys.path.append("../FidlTrack_example_data/240130_cos418+716_3.5ul_6ms")

from config_tracking import *

filenames = []
for root, dirs, files, in os.walk(base_dir):
	if not any([e in root for e in exclude]):
		filenames.extend(["/".join([root, f]) for f in files if f.endswith(base_ext) and f.startswith(base_start) and all([e not in f for e in exclude])])


for cpt, fname in enumerate(filenames):
	print("Processing[{}/{}]: {}".format(cpt+1, len(filenames), fname))
	exp_path = "/".join(fname.split("/")[:-1])
	exp_fname = fname.split("/")[-1]

	base_fname = exp_fname[len(base_start):-len(base_ext)]

	if mask:
		mask_f = "/".join([exp_path, mask_fname.format(fname=base_fname)])
		print(" Masking spots: " + mask_f)
		if not path.isfile(mask_f):
			print(" [ERROR] Mask file not found")
			continue

	outPath = "/".join([out_dir, base_fname])

	if not path.isdir(outPath):
		os.makedirs(outPath)

	IJ.open(fname)
	imp = IJ.getImage()
	imp_name = imp.getTitle()
	dims = imp.getDimensions()

	mask_imp = None
	if mask:
		IJ.open(mask_f)
		mask_imp = IJ.getImage()
		IJ.selectWindow(imp_name)
		imp = IJ.getImage()
		is_single_frame = all([e == 1 for e in mask_imp.getDimensions()[2:]])

	dims = imp.getDimensions()
	if dims[3] > 1:
		print("Corrected T and slice dimensions")
		imp.setDimensions(dims[2], dims[4], dims[3])

	if do_calibration:
		print("Performing calibration")
		calib = Calibration()
		calib.setTimeUnit("sec")
		calib.setUnit("micron")
		if set_dt is not None:
			calib.frameInterval = set_dt
		if set_dx is not None:
			calib.pixelWidth = set_dx
			calib.pixelHeight = set_dx
		imp.setCalibration(calib)

	calib = imp.getCalibration()
	print(" Pixel size: {} um, Acquisition time: {} s".format(
		calib.pixelWidth, calib.frameInterval))

	#===== SPOT DETECTION
	for p_DIAMETER in p_DETECTION_BLOB_DIAMETERS:
		cur_title = imp.getTitle()
		imp = WindowManager.getImage(cur_title)

		for spot_th in ths:
			out_fname = "/".join([outPath, spot_fname.format(mask=mask, spot_rad=p_DIAMETER, spot_th=spot_th)])
			if not force and path.isfile(out_fname):
				print("  Skipped")
				continue
			print("  {}".format(out_fname))

			imp.show()

			model = Model()
			model.setLogger(Logger.IJ_LOGGER)

			settings = Settings(imp)

			# Configure spot detector
			settings.detectorFactory = LogDetectorFactory()
			settings.detectorSettings = { 
			    'DO_SUBPIXEL_LOCALIZATION': p_DO_SUBPIXEL_LOCALIZATION,
			    'RADIUS': p_DIAMETER / 2, #DIAMETER TO RADIUS .....................
			    'THRESHOLD': float(spot_th),
			    'DO_MEDIAN_FILTERING': p_DO_MEDIAN_FILTERING,
			    'TARGET_CHANNEL': 1
			} 

			#remove spots too close to the field of view boundary that cause artifacts
			settings.addSpotFilter(FeatureFilter('POSITION_X', p_DIAMETER / 2, True))
			settings.addSpotFilter(FeatureFilter('POSITION_X', dims[0] * calib.pixelWidth - p_DIAMETER / 2, False))
			settings.addSpotFilter(FeatureFilter('POSITION_Y', p_DIAMETER / 2, True))
			settings.addSpotFilter(FeatureFilter('POSITION_Y', dims[1] *  calib.pixelWidth - p_DIAMETER / 2, False))

			trackmate = TrackMate(model, settings)
			trackmate.setNumThreads(4)

			ok = trackmate.execDetection()
			if not ok:
				print(str(trackmate.getErrorMessage()))
				continue
			ok = trackmate.computeSpotFeatures(True)
			if not ok:
				print(str(trackmate.getErrorMessage()))
				continue
			ok = trackmate.execSpotFiltering(True)
			if not ok:
				print(str(trackmate.getErrorMessage()))
				continue

			if "reg" in locals() and reg:
				print(" REGISTERING SPOTS: {},{}".format(reg[0], reg[1]))

			if mask_imp:
				n = 0
				spots_to_rm = []
				for s in model.getSpots().iterator(True):
					n += 1
					pos = [s.getFeature("POSITION_X"), s.getFeature("POSITION_Y")]

					if "reg" in locals() and reg:
						pos = [pos[0] + reg[0], pos[1] + reg[1]]
						s.putFeature("POSITION_X", pos[0])
						s.putFeature("POSITION_Y", pos[1])

					p = [int(s.getFeature("FRAME")), int(pos[0] // set_dx), int(pos[1] // set_dx)]
					if is_single_frame:
						if mask_imp.getStack().getProcessor(1).getPixel(p[1], p[2]) == 0:
							spots_to_rm.append(s)
					elif mask_imp.getStack().getProcessor(p[0]+1).getPixel(p[1], p[2]) == 0:
						spots_to_rm.append(s)

				base_spts = model.getSpots()
				print(" REMOVING {}/{} spots not in structure".format(len(spots_to_rm), n))
				for s in spots_to_rm:
					model.removeSpot(s)
			
			###REMOVE DUPLICATED SPOTS (WHY IS THIS EVEN A THING???)
			spts = list(model.getSpots().iterator(True))
			spts_fr = {}
			for i in range(len(spts)):
				fr = int(spts[i].getFeature("FRAME"))
				if fr not in spts_fr:
					spts_fr[fr] = []
				spts_fr[fr].append(i)

			to_rm = []
			for i in range(len(spts)):
				f1 = int(spts[i].getFeature("FRAME"))
				for j in spts_fr[f1]:
					if i == j:
						continue
					f2 = int(spts[i].getFeature("FRAME"))
					if f1 != f2:
						continue
					if sqrt(spts[i].squareDistanceTo(spts[j])) < 1e-7:
						if spts[i].getFeature("QUALITY") > spts[j].getFeature("QUALITY"):
							to_rm.append(spts[j])
						else:
							to_rm.append(spts[i])
			print(" REMOVING {}/{} duplicated spots".format(len(to_rm), len(spts)))
			for s in to_rm:
				model.removeSpot(s)

			towrite = ["POSITION_T", "POSITION_X", "POSITION_Y", "FRAME", "RADIUS", "QUALITY"]
			with open(out_fname, 'w') as f:
				f.write(",".join(towrite) + "\n")
				for spot in model.getSpots().iterator(True):
					feats = spot.getFeatures()
					f.write(",".join([str(feats[k]) for k in towrite]) + "\n")

			model.clearSpots(True)

	if mask_imp:
		mask_imp.changes = False
		mask_imp.close()

	imp.changes = False
	imp.close()

print("done")
