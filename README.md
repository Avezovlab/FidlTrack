# FidlTrack, a platform for improved single-particle tracking fidelity.

FidlTrack is composed of three modules:

## Module 1: Experimental and tracking setup optimiser

Finds the fidelity-maximising spot density and linking distance
parameters, given an estimate of the diffusion coefficient of the
particle of interest and the acquisition framerate. This module works
either in freespace (no structural constraint on motion, e.g. plasma
membrane), or in ER or mitochondria network geometries.

The optimiser is provided as a Python notebook `notebooks/FidlTrack_predict.ipynb`, also accessible through Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Avezovlab/FidlTrack/blob/main/notebooks/FidlTrack_predict.ipynb)


## Module 2: Structure-aware tracking

Inserts geometrical constraint into tracking algorithms to improve
tracking accuracy. The pipeline to use structure-aware tracking is as
follow:

0. Before starting, you need to have already performed spot detection
on your data, this can be achieved using the imageJ script
imageJ_tracking_scripts/trackmate_spot_detection_multi_th py.

1. obtain a mask of the structure constraining the particle motion:
either pre/post single-particle acquisition or simultaneously in 2
colours. The structure-aware works both with a single static mask or a
stack of mask.

2. Process the mask stack (or image) using the Python notebook
`colab_notebooks/FidlTrack_buildStructureGraph.ipynb`, also accessible
through Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Avezovlab/FidlTrack/blob/main/notebooks/FidlTrack_buildStructureGraph.ipynb)

This notebook transform the mask(s) into a pre-computed distance file
along each of the provided structures.

3. Fill the config_tracking.py file, a template is provided in this
folder. You will need to retrieve the following files from step 2:
    + windowed component stack (or component image for single mask
    input) typically ending with "_comps_wdur=XX_wovlp=YY.tif"
    + the pre-computed graph distances stored in an "optimised" binary
    file format typically ending with "_dist=XX.bin"

Run the python script
imageJ_tracking_scripts/trackmate_struct_aware_tracking.py inside
imageJ using the config_tracking.py file as input.

4. The resulting tracking files are similar to conventional trackmate
track files with extra column providing structure-aware information:
      * winIdx: time-window index of the spot.
      * component: connected component associated to the spot.
      * d_g (um) : graph distance.
      * Ambiguity : number of possible successor spots for this displacement
    (see module 3).


## Module 3: Ambiguity detection, quantification and removal

A trajectory displacement, a segment connecting two succesive points
A,B inside a trajectory, is ambiguous when there are more than one
possible valid spot to choose from for point B. Technically, it means
there are multiple spots whithin the linking distance of A in its next
frame. Ambiguous displacements are particularly error-prone as the
tracking algorithm had to make a choice that we have no guarantee of
being locally optimal. The Ambiguity Score is the percentage of
ambiguous displacements found in all trajectories and is a good
measure of the fidelity of a dataset, as according to our simulations
> 50% of tracking errors happen at ambiguous displacements. Note
however that ambiguity is a conservative measure and not all ambiguous
displacements are necessarily erroneous.

The following codes are provided to handle ambiguities:

1. We developed an extension to the Trackmate imageJ plugin that
provides an EdgeAmbiguityAnalyzer that adds an extra "Ambiguity"
column to tracking files. This column reports the number of possible
successor for each displacement minus one (so that 0 is a
non-ambiguous displacement and > 0 is ambiguous), the last
displacement of each trajectory is given the value -2 (no ambiguity
can be computed here as there is no successor).

You can run it for conventional tracking using the imageJ script
provided in imageJ_tracking_scripts/trackmate_tracking_from_spots.py.

The EdgeAmbiguityAnalyzer is also directly integrated in the
structure-aware tracking extension and where it reports ambiguities
based on the graph distance (instead of Euclidean distance for the
conventional tracking).


2. Once ambiguities are detected, we provide a Python notebook
`notebooks/FidlTrack_ambiguity.ipynb`, to quantify and remove ambiguous
displacements displacements, also accessible through Google Colab:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Avezovlab/FidlTrack/blob/main/notebooks/FidlTrack_ambiguity.ipynb)

This notebook is decomposed in three parts:

* Inputs, where you can provide either :
    * a tracking files with an Ambiguity column (given by the
    EdgeAmbiguityAnalyzer);
    * a tracking file without ambiguity information AND the corresponding detected spot
    file, in which case ambiguous displacements are computed (this option takes longer).

* Then, the Ambiguity Score is computed giving the percentage of
  ambiguous displacements in the trajectories. This score can be used
  as a quality measure of a dataset.

* Finally, the notebook provides the possibility to remove ambiguous
  displacements, by splitting affected trajectories into two
  sub-trajectories before and after the displacements. This procedure
  improves the reliability of a dataset at the cost of fragmenting
  trajectories (increases the amount of trajectories and decreases
  their lengths).

If working inside complex structures (e.g. neurites, ER,
mitochondria), ambiguity is best coupled with structure-aware tracking
as the graph distance strongly decreases the amount of ambiguous
displacements.



