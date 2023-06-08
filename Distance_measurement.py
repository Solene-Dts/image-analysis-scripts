###Introduction

# Hi!

# Characteristics
	# written for Jython (Python 2.7 in ImageJ 2.9.0)
	# contact for help: solenedeutsch@gmail.com

# Aim of the script
	# Isolates vessels and T cells and measures the distance from each T cell to each vessel.

# Getting started
	# open the images you want to analyze with the "Split channel" box unchecked
		# the results of one run will all be put in the same text file
	# the results are written in a text file, in your choosen directory
	# the results are formatted to be compatible with the distance_analysis script to analyze them
	# a preliminary run with fewer images is required to make sure that all the thresholdings and detections are working correctly
	# this script selects the active image, analyzes it, and closes it. Make sure no useless image is opened in front.
	# the T cell detection works with Cellpose. See Cellpose-and-Cellpose-on-ImageJ Tutorial for more information.
	# you can change the name of the file at the bottom of the script (Saving results part)
	# This script currently asks user to check that the T cells detected are indeed T cells (and same for the vessels).
		# If you wnat to run the script without having to check, comment the lines where you see a WaitForUserDialog (2 of them)
	
# Have fun!

###Imports

from ij import IJ, Prefs
from ij.plugin.frame import RoiManager
from ij.plugin import ChannelSplitter, ZProjector, RoiEnlarger
from ij.gui import WaitForUserDialog
from ij.measure import ResultsTable

###Preliminary functions

def define_ROI(vessels, OriImage, show = False):
	'''
	Takes the z-projection and returns a list of z-stacks, each containing a single vessel.
	Both parameters are images, returns a list of images.
	'''
	# Isolate individual vessels
	rm = RoiManager.getRoiManager()
	n_roi = rm.getCount() # number of ROIs (=number of vessels detected)
	list_vessels = [] 
	ra = rm.getRoisAsArray()
	# Apply ROI and clear outside to keep only the vessel part
	for i in range(n_roi):
		single_vessel = OriImage.duplicate()
		single_vessel.setRoi(ra[i])
		IJ.setBackgroundColor(0,0,0)
		IJ.run(single_vessel, "Clear Outside", "stack")
		IJ.run(single_vessel, "Select None", "stack")
		list_vessels.append(single_vessel)
		if show:
			single_vessel.show()
			check = WaitForUserDialog("Time to look at the data", "Check vessel.")
			check.show()
			single_vessel.changes = False
			single_vessel.close()
	return list_vessels

def retrieve_data(isolate_vessels, marker_image, image_name):
	'''
	Retrieves information about the size of the vessel, and the intensity of the marker in marker_image in that vessel.
	'''
	ra = rm.getRoisAsArray()
	# Isolate the shape of the vessel in the image containing the marker staining
	list_intensities = define_ROI(isolate_vessels, marker_image)
	image_to_analyze = list_intensities[int(image_name[-1])]
	# Z-projection
	image_to_analyze = ZProjector.run(image_to_analyze,"max")
	# Setting the ROI around the vessel
	image_to_analyze.setRoi(ra[int(image_name[-1])])
	# If it has't already been done, retrieves the size of the vessel
	if len(DATA[image_name]) == 0:
		DATA[image_name].append(image_to_analyze.getStatistics().area)
	# Retrieves the intensity of the staining inside the vessel
	DATA[image_name].append(image_to_analyze.getStatistics().mean)

def finding_Tcells(t_chan):
	'''
	Takes the image containing CD3 staining and isolates the T-cells.
	Stores Rois containing the T-cells in the Roi Manager and returns the number of T-cells found.
	'''
	t_cells = t_chan.duplicate()
	# Image treatment
	IJ.setAutoThreshold(t_cells, "Moments dark")
	Prefs.blackBackground = True
	IJ.run(t_cells, "Convert to Mask", "method=Moments background=Dark calculate black")
	t_cells = ZProjector.run(t_cells,"sum")
	t_cells.show()
	rm = RoiManager.getRoiManager()
	not_empty = len(rm.getRoisAsArray())
	if not_empty != 0:
		rm.runCommand(t_cells, "Delete")
	# Running cellpose for T cell detection
	IJ.run("Cellpose Advanced (custom model)", "diameter=0 cellproba_threshold=-1.0 flow_threshold=0.6 anisotropy=1.0 diam_threshold=12.0 model_path=C:\Users\solde212\.cellpose\models\Tcells_GL261 model=Tcells_GL261 nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=0 omni=false cluster=false additional_flags=")
	# Retrieving proper ROIs from cellpose
	# cellpose outputs a mask that contains found areas, each filled with a different pixel intensity value (starts at 1, step=1)
	mask = IJ.getImage()
	IJ.run(mask, "Select All", "")
	roi = mask.getRoi()
	rm.addRoi(roi)
	table = ResultsTable.getResultsTable()
	IJ.run("Clear Results", "")
	rm.runCommand(mask,"Measure")
	max_int = table.getValue("Max", 0)
	total_roi = range(1, int(max_int)+1)
	rm.setSelectedIndexes([0])
	rm.runCommand(mask, "Delete")
	# Isolating T-cells
	for value in total_roi:
		temp = mask.duplicate()
		IJ.setThreshold(temp, value, value)
		IJ.run(temp, "Convert to Mask", "")
		IJ.run(temp, "Create Selection", "")
		roi = temp.getRoi()
		rm.addRoi(roi)
		temp.close()
	# Clearing useless windows for the next image
	IJ.selectWindow("Results")
	IJ.run("Close")
	mask.close()
	# Adding Roi names for clarity (useful for debugging, can be commented out otherwise)
	for i in xrange(len(rm.getRoisAsArray())):
		rm.select(i)
		rm.runCommand("Rename", "T_"+str(i))
	# Manually checking that the ROIs are actuall T-cells
	check = WaitForUserDialog("Time to look at the data", "Do those look like T cells ? \nIf not, modify the ROIs then click OK.")
	check.show()
	t_cells.close()
	return len(rm.getRoisAsArray())

def distance_to_vessels(vessel, vessel_key):
	'''
	Performs the distance analysis.
	vessel_stack : image containing the vessel
	vessel_key : code name of the vessel image
	'''
	# Analysis is performed on a z projection
	rm = RoiManager.getRoiManager()
	nb_tcells = len(rm.getRoisAsArray())
	# Selecting the vessel inner area
	IJ.setAutoThreshold(vessel, "Li dark")
	Prefs.blackBackground = True
	IJ.run(vessel, "Convert to Mask", "")
	IJ.run(vessel, "Create Selection", "")
	# Adding vessel inner area Roi to the analysis
	roi = vessel.getRoi()
	rm.addRoi(roi)
	# Adding Roi names for clarity (useful for debugging, can be commented out otherwise)
	rm.select(nb_tcells)
	rm.runCommand("Rename", "Area_0")
	# Initializing parameters for the loop
	previous = nb_tcells
	current = nb_tcells +1
	new = nb_tcells +2
	# Creating distance layers
	for i in xrange(5):
		# Roi containing the vessel, expended by 10µm each time
		rm.select(previous)	
		roi = RoiEnlarger.enlarge(roi, 55)
		rm.addRoi(roi)
		rm.runCommand(vessel,"Deselect")
		rm.setSelectedIndexes([previous, current])
		rm.runCommand(vessel,"XOR")
		new_roi = vessel.getRoi()
		rm.addRoi(new_roi)
		rm.runCommand(vessel,"Deselect")
		rm.select(new)
		rm.runCommand("Rename", "Area_"+str(i+1))
		# Updating loop parameters
		previous = current
		current +=2
		new += 2
	# Add Roi containing evrything except the largest layer's area
	IJ.run(vessel, "Select All", "")
	roi = vessel.getRoi()
	rm.addRoi(roi)
	rm.setSelectedIndexes([previous, current])
	rm.runCommand(vessel,"XOR")
	new_roi = vessel.getRoi()
	rm.addRoi(new_roi)
	rm.runCommand(vessel,"Deselect")
	rm.select(new)
	rm.runCommand("Rename", "Area_"+str(i+2))
	# Getting the position of each T-cell found
	for tcell in xrange(nb_tcells):
		list_overlap_areas = []
		# Checking if the T-cell is in any of the layers previously defined
		for area in xrange(nb_tcells, len(rm.getRoisAsArray()), 2):
			rm.setSelectedIndexes([tcell,area])
			rm.runCommand(vessel,"AND") # Checking if a Roi that would contain both a part of the T-cell and a part of the layer exists
			if vessel.getRoi() != None:
				# If it exists, calculates its area
				new_roi = vessel.getRoi()
				rm.addRoi(new_roi)
				roiStat = new_roi.getStatistics()
				overlap_area = roiStat.area
				list_overlap_areas.append((area-nb_tcells)/2)
				list_overlap_areas.append(overlap_area)
		# If only one overlap found between the T-cell and a layer...
		if len(list_overlap_areas) ==2:
			#...keep the number of the layer
			DATA[vessel_key].append(list_overlap_areas[0])
		# If a T-cell overlaps with more than 1 layer...
		elif len(list_overlap_areas) >2:
			#...keep the layer that contains most of the T-cell
			size1 = list_overlap_areas[1]
			size2 = list_overlap_areas[3]
			if size1 > size2:
				DATA[vessel_key].append(list_overlap_areas[0])
			else:
				DATA[vessel_key].append(list_overlap_areas[2])
	for value in range(len(DATA[vessel_key])):
		# A distance value of 6 means the T cell was outside all layers
		if DATA[vessel_key][value] == 6:
			DATA[vessel_key][value] = "too far"
	# Manually checking if the ROI is actually a vessel
	myWait = WaitForUserDialog("Time to look at the data", "Does this region look like a vessel ? \nIf not, write it down then click OK.\nThe data concerning the fake vessel can then be removed by hand")
	myWait.show()
	# Close useless windows and ROIs
	ra = rm.getRoisAsArray()
	for i in reversed(range(n_tcells, len(ra))):
	 	rm.setSelectedIndexes([i])
	 	rm.runCommand(whole, "Delete")
	vessel.changes = False
	vessel.close()
	
###Initializing

rm = RoiManager.getRoiManager()

DATA = {}

count = int(IJ.getNumber("How many images to analyze ?", 1))

Dir = os.getcwd()
intermediate = IJ.getString("Write directory for the result file (leave empty if current directory is wanted) :", "E:\\PROJECTS\Solene")
if intermediate == "":
	Wdir = Dir
else :
	Wdir = intermediate

vessel_channel_number = int(IJ.getString("Enter number of channel for vessel staining", "2"))
marker1_channel_number = int(IJ.getString("Enter the number of the channel you are interested in measuring", "1"))
tcell_channel_number = int(IJ.getString("Enter the number of the T-cell channel", "0"))

tumor_number = IJ.getString("Enter tumor number (for file naming)", "")
	
### Analysis

for _ in xrange(count):
	# Initializing for the current image
	
	intensity = {}
	
	whole = IJ.getImage()
	
	channels = ChannelSplitter.split(whole)
	
	tcell_channel = channels[tcell_channel_number]
	vwf_channel = channels[marker1_channel_number]
	vessel_channel = channels[vessel_channel_number]
	
	# Choosing the vessels
	
	vessels = ZProjector.run(vessel_channel,"max")
	ori_vessel = vessels.duplicate()
	IJ.setAutoThreshold(vessels, "Moments dark")
	Prefs.blackBackground = True
	IJ.run(vessels, "Convert to Mask", "")
	IJ.run(vessels, "Despeckle", "")
	IJ.run(vessels, "Despeckle", "")
	IJ.run(vessels, "Despeckle", "")
	IJ.run(vessels, "Options...", "iterations=1 count=1 black do=Dilate")
	IJ.run(vessels, "Options...", "iterations=7 count=1 black do=Close")
	IJ.run(vessels, "Analyze Particles...", "size=30-Infinity show=Nothing clear add")
	show = False
	#show = True #Decomment to follow vessel detection
	list_vessels = define_ROI(vessels, ori_vessel, show)
	#Computing data for each vessel in the image
	for k in xrange(len(list_vessels)):
		vessel = list_vessels[k]
		v_title = tumor_number+"_"+vessel.getTitle()[-1]+"_"+str(k)
		DATA[v_title] = []
		#Get the size and intensity of the staining
		retrieve_data(vessel, vwf_channel, v_title)
	#Detect T cells in the image
	n_tcells = finding_Tcells(tcell_channel)
	#Calculate the distance
	for s in xrange(len(list_vessels)):
		vessel = list_vessels[s]
		v_title = tumor_number+"_"+vessel.getTitle()[-1]+"_"+str(s)
		distance_to_vessels(vessel, v_title)
	
	# Closing current image
	
	whole.close()

#Quick printing for debugging
print(DATA)
	
### Saving the data

results = open(Wdir+"/new_distance_"+tumor_number+".txt", "w") 
results.write(str(DATA))
results.close()		
	