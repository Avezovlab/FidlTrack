force = False #this will force recomputing all results

#calibrate the files (set pixelsize and acquisition times)
do_calibration = False
set_dx = XX
set_dt = XX

#directories
base_dir = XX
out_dir = XX

#These three fields allow to find the files to process
base_start = XX #start of the file, use "" to skip
base_ext = XX #extension of the file, use "" to skip
exclude = [] #list of names to exclude


###The spot detection and tracking scripts will process all the files in a sub-folder of base_dir
###starting with base_start and ending with base_ext
###The file being processed, without start and extension (fname[len(base_start):-len(base_ext)],
###will be passed to mask_fname, comp_fname_win and dist_fname as the variable fname

#struct stack mask filtering
mask = False
is2D = False
mask_fname = "XX{fname}XX" #compute the mask filename from the spot filename
comp_fname_win = "XX{fname}XX" #compute the component filename from the spot filename
dist_fname = "XX{fname}XX" #compute the distance filename from the spot filename

##Spot detection
#Trackmate
p_DO_SUBPIXEL_LOCALIZATION = True
p_DO_MEDIAN_FILTERING = True
p_DETECTION_BLOB_DIAMETERS = [XX] #list of spot diameter values to consider
ths = [XX] #list of threshold values to consider
spot_fname = "spots_mask={mask}_rad={spot_rad}_th={spot_th}.csv"
reg = [] #registration x,y coordinates, leave empty or None if not needed

###For the following variables, fname is the name of the spot file without csv extension

##Tracking
p_LINKING_MAX_DISTANCES = [XX] #list of linking distance values to consider
p_MAX_FRAME_GAPS = [0] #0 for structure-aware tracking
link_fname = "tracks_{fname}_dist={dist}_distgap={distgap}_framegap={framegap}.csv"  #for conventional tracking
link_struct_fname = "tracks_struct_{fname}_dist={dist}_distgap={distgap}_framegap={framegap}.csv" #for struct-aware tracking
