###Introduction

# Hi!

# Characteristics
	# written for Jython (Python 2.7 in ImageJ 2.9.0)
	# contact for help: solenedeutsch@gmail.com

# Aim of the script
	# Isolates vessels and measures their area and the corresponding intensity of two different markers.

# Getting started
	# open the images you want to analyze with the "Split channel" box unchecked
		# /!\ open only what you want to analyze, the results of one run will all be put in the same text file
	# the results are written in a text file, in your choosen directory
	# the results are formatted to be compatible with the subcluster_graphs script to analyze them
	# a preliminary run with fewer images is required to make sure that all the thresholdings and detections are working correctly
	# this script selects the active image, analyzes it, and closes it. Make sure no useless image is open.
	
# Have fun!

###Imports

from ij import IJ, Prefs
from ij.plugin.frame import RoiManager
from ij.plugin import ChannelSplitter, ZProjector
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

def retrieve_data(isolate_vessels, marker_image, image_name, threshold_min):
	'''
	'''
	list_intensities = define_ROI(isolate_vessels, marker_image)
	image_to_analyze = list_intensities[int(image_name[-1])]
	rm.setSelectedIndexes([int(image_name[-1])])
	rm.runCommand(image_to_analyze,"Measure")
	IJ.setThreshold(image_to_analyze, threshold_min, 255)
	Prefs.blackBackground = True
	IJ.run(image_to_analyze, "Convert to Mask", "")
	IJ.run(image_to_analyze, "Divide...", "value=255.000")
	rm.setSelectedIndexes([int(image_name[-1])])
	rm.runCommand(image_to_analyze,"Measure")
	
def set_threshold(image):
	IJ.run(image, "Select All", "")
	roi = image.getRoi()
	rm.addRoi(roi)
	ra = rm.getRoisAsArray()
	rm.setSelectedIndexes(range(len(ra)))
	rm.runCommand(image,"XOR")
	threshold = image.getStatistics().mean
	rm.runCommand(image,"Deselect")
	rm.setSelectedIndexes([len(ra)])
	rm.runCommand(image, "Delete")
	return threshold

	
###Initializing

rm = RoiManager.getRoiManager()

DATA = {}

count = int(IJ.getNumber("How many images to analyze ?", 1))

Dir = "E:\\PROJECTS\Solene"
intermediate = IJ.getString("Write directory for the result file (leave empty if current directory is wanted) :", "E:\\PROJECTS\Solene\P1SD12")
if intermediate == "":
	Wdir = Dir
else :
	Wdir = intermediate

vessel_channel_number = int(IJ.getString("Enter number of channel for vessel staining", "3"))
marker1_channel_number = int(IJ.getString("Enter the number of the channel your interested in measuring", "1"))
marker2_channel_number = int(IJ.getString("Enter the number of the other channel your interested in", "2"))

experiment_type = IJ.getString("Enter experiment type (will be part of the result file's title)", "tumor")
tumor_number = IJ.getString("Enter tumor number", "")
	
### Analysis

for _ in xrange(count):
	# Initializaing for the current image
	
	intensity = {}
	
	whole = IJ.getImage()
	
	channels = ChannelSplitter.split(whole)
	
	plvap_channel = channels[marker1_channel_number]
	vwf_channel = channels[marker2_channel_number]
	vessel_channel = channels[vessel_channel_number]
	
	# Choosing the vessels
	
	vessels = ZProjector.run(vessel_channel,"max")
	vwf = ZProjector.run(vwf_channel,"max")
	plvap = ZProjector.run(plvap_channel,"max")
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
	threshold_plvap = set_threshold(plvap)
	threshold_vwf = set_threshold(vwf)
	for k in xrange(len(list_vessels)):
		vessel = list_vessels[k]
		v_title = vessel.getTitle()[-12:]+"_"+str(k)
		DATA[v_title] = []
		retrieve_data(vessel, plvap, v_title, threshold_plvap)
		retrieve_data(vessel, vwf, v_title, threshold_vwf)
		table = ResultsTable.getResultsTable()
		vessel_area = [area for area in table.getColumn(0)][0]
		vessel_intensities = [intensity for intensity in table.getColumn(1)]
		DATA[v_title].append(vessel_area)
		DATA[v_title] += vessel_intensities
		IJ.run("Clear Results", "")
	
	# Closing current image
	
	whole.close()

print(DATA)
	
# Saving the data

results = open(Wdir+"/"+experiment_type+"_"+tumor_number+".txt", "w") 
results.write(str(DATA))
results.close()		
	
	
	
	
	