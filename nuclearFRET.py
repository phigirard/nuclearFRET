#*******************************************************************************
#                                                                             
#	Philippe GIRARD 
# 	Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
# 	nuclearFRET.py
#	Release v1.0
#
#	Copyright 2022 - AGPL-3.0 License
#                                                                             
#******************************************************************************/


#@ File(label="LSM Multidimensional file", description="Select a Hyperstack CTZ-(Series or not) file ", style="file") inputFile
#@ String (label="Type of FRET analysis", description="Select the object to analyze ", choices={"Nuclei", "Whole cell"}, style="radioButtonHorizontal") FRETtype
#@ Boolean(label="Choose values for background subtraction", description="Values for background subtraction",value=False, persist=True) backgroundSubtract
#@ Double(label="Background intensity level of Donor channel", description="Background intensity level of donor?",value=140, persist=True)  backgroundDonor
#@ Double(label="Background intensity level of Acceptor channel", description="Background intensity level of acceptor?",value=140, persist=True)  backgroundAcceptor



#@ UIService uiService
#@ LogService log

#---------------------------------------------------------------
#---------------          LIBRARIES            -----------------
#---------------------------------------------------------------

# Python Library -------------------------------------------------------------------------------------------
import os
from os import path
from collections import OrderedDict

import sys


# ImageJ Library ----------------------------------------------------------------------------
from ij import IJ, ImagePlus,  Prefs, WindowManager
from ij.process import ImageConverter, ByteProcessor, AutoThresholder,  FloatPolygon, ImageStatistics
from ij.gui import GenericDialog, WaitForUserDialog, PlotWindow, ProfilePlot, Overlay, Line, Wand
from ij.gui import Plot, Roi, PointRoi,PolygonRoi, OvalRoi, ShapeRoi
from ij.plugin import ZAxisProfiler, ZProjector
from ij.plugin import RoiEnlarger
from ij.plugin import ImageCalculator as IC
from ij.measure import ResultsTable , Measurements, Calibration
from ij.plugin.filter import Analyzer
from ij.plugin.frame import RoiManager
from ij.text import TextWindow

# Loci Library ------------------------------------------------------------------------------------------
from loci.formats import ImageReader
from loci.formats import MetadataTools
from loci.plugins.in import ImporterOptions
from loci.plugins import BF


# MorphoLibJ Library ---------------------------------------------------------------------------------------
from inra.ijpb.watershed import Watershed, MarkerControlledWatershedTransform2D
from inra.ijpb.binary import BinaryImages
from inra.ijpb.label import LabelImages
from inra.ijpb.color.ColorMaps import CommonLabelMaps
from inra.ijpb.color import CommonColors

# Java class ---------------------------------------------------------------------------------------
from java.lang import Float
from java.awt import Color
from java.awt.event import AdjustmentListener  
from java.awt.image import BufferedImage
from java.awt.geom import Rectangle2D, Ellipse2D

from  javafx.geometry import Point2D
from  math import isnan

from org.jfree.chart import ChartPanel, JFreeChart
from org.jfree.chart.axis import NumberAxis
from org.jfree.chart.plot import XYPlot
from org.jfree.chart.renderer.xy import XYLineAndShapeRenderer
from org.jfree.data.xy import XYDataset, XYSeries, XYSeriesCollection

#---------------------------------------------------------------
#---------------          CONSTANTS            -----------------
#---------------------------------------------------------------
setconcat = False
showomexml = False
autoscale = False
Prefs.blackBackground = True


showRoiManager = False # show RoiManager

FRETTitle = "FRET index (%)" #Title of the RT for FRET measurements

IJ.setForegroundColor(255, 255, 255) # set foreground color to white


#---------------------------------------------------------------
#----------------- All Classes for analysis  -----------------
#---------------------------------------------------------------


## threshold class
class ThresholdPreviewer(AdjustmentListener):  
	def __init__(self, imp, sliders):  
		self.imp = imp   
		#self.currentSlice = imp.getCurrentSlice()
		self.sliders = sliders  
		self.threshold()
		
	def adjustmentValueChanged(self, event):  
		self.threshold()   
			
	def reset(self):  
		IJ.resetThreshold(self.imp) 
		
	def threshold(self):
		IJ.setThreshold(self.imp, self.sliders.get(0).getValue() ,self.sliders.get(1).getValue())
	
	def getMinThreshold(self) :
		return self.sliders.get(0).getValue()
		
	def getMaxThreshold(self) :
		return self.sliders.get(1).getValue()

#---------------------------------------------------------------
#----------------- All Functions for analysis  -----------------
#---------------------------------------------------------------


#### Fonctions for STEP 1 : Preparation of data & analysis - Preprocessing

# Get indexes of Donnor and Acceptor images
def getImpIndexes(imagefile_ , sizeC_ , idxSeries_ ):
	#read in ImagePlus with arguments
	options = ImporterOptions()
	options.setId(imagefile_)
	options.setSeriesOn(idxSeries_,True)	
	impsLSM = BF.openImagePlus(options)

	#projection in Z dimension :impProj_
	IJ.run(impsLSM[0], "Grays", "stack")
	impProj_ = ZProjector.run(impsLSM[0],"max")
	IJ.run(impProj_, "Enhance Contrast", "saturated=0.35")	
	IJ.setTool("rectangle")
	impProj_.show()
	idProj_ = impProj_.getID()
	IJ.run("Brightness/Contrast...")
	waitDialog = WaitForUserDialog("ROI","Select a ROI with the rectangle tool")
	waitDialog.show()

	impsLSM[0].setCalibration(Calibration())
	impsLSM[0].setDimensions(1, 1, sizeC_)
	plot = ZAxisProfiler.getPlot(impsLSM[0])
	plot.setXYLabels('Channel','Mean')
	xvalues = plot.getXValues()
	yvalues = plot.getYValues()
	series = XYSeries("Mean intensity channel")
	for i in range(len(xvalues)) :
		series.add(xvalues[i], yvalues[i])
	dataset = XYSeriesCollection(series) #Add series to dataset
	yaxis = NumberAxis("Mean")
	xaxis = NumberAxis("# Channel")
	offset = 0.5
	xaxis.setRange(1-offset, sizeC_+offset)
	r = XYLineAndShapeRenderer()
	r.setSeriesPaint(0, Color.RED)
	r.setSeriesShape(0, Ellipse2D.Double(-3.0,-3.0,6.0,6.0))
	xyplot = XYPlot(dataset, xaxis, yaxis, r)
	xyplot.setBackgroundPaint(Color.white)
	chart = JFreeChart(xyplot)
	chart.removeLegend()
	impPlot = IJ.createImage("Mean intensity channel", "RGB", 512, 512, 1);
	imagePlot = impPlot.getBufferedImage()
	chart.draw(imagePlot.createGraphics(), Rectangle2D.Float(0, 0, impPlot.getWidth(), impPlot.getHeight()))
	impPlot.setImage(imagePlot)
	impPlot.show()

	gui = GenericDialog("Select FRET Donnor/Acceptor images")
	gui.addSlider("Donnor image: ", 1, sizeC_, 3)
	gui.addSlider("Acceptor image: ", 1, sizeC_, 6)
	gui.showDialog()
	impPlot.close()
	impProj_.deleteRoi()
	idxDonnor_ = int(gui.getNextNumber())
	idxAcceptor_ = int(gui.getNextNumber())
	for i in range(len(impsLSM)) :
		impsLSM[i].flush()
		
	return idxDonnor_, idxAcceptor_, idProj_

	
#Open ZCT hyperstack with only one channel index 
def extractImpFromIndex(imagefile_ , idxChannel_ , idxSeries_ ):
	options = ImporterOptions()
	options.setId(imagefile_)
	options.setSeriesOn(idxSeries_,True)	
	options.setAutoscale(autoscale)
	options.setShowOMEXML(showomexml)
	options.setConcatenate(setconcat)
	options.setCBegin(idxSeries_ , idxChannel_ -1)
	options.setCEnd(idxSeries_ ,idxChannel_ -1)  
	return BF.openImagePlus(options)[0]

# Adjust the string size with 0
def adjustSizeNum(s_ , length_):
	while len(s_)<length_ :
		s_ = "0"+s_
	return s_

#create folder for each series
def createFolder(imagefile_ , idxSeries_):
	srcDir_ = os.path.dirname(imagefile_)
	basename_ = os.path.basename(os.path.splitext(imagefile_)[0]).replace(' ', '_').lower()
	basename_ += "_S" + adjustSizeNum(str(idxSeries_), 2)
	anDir_ = os.path.join(srcDir_, basename_) 
	if not os.path.exists(anDir_):
		os.makedirs(anDir_)
	return anDir_

# select the background ROI from an image 
def getBackgroundROI(imp_) :
	print 'Select area in the background'
	imp_.show()
	IJ.run(imp_, "Enhance Contrast", "saturated=0.35")
	IJ.setTool("freeline")
	myWait = WaitForUserDialog("Background selection", "Select a ROI in the background with the freeline tool")
	myWait.show()
	roi_  = imp_.getRoi()
	imp_.deleteRoi()
	imp_.hide()
	return roi_

# subtract the mean intensity of the ROI obtained from the function getBackgroundROI
def subtractBackground(imp_, roi_):
	Analyzer.setMeasurements(Measurements.MEAN)
	Analyzer.setPrecision(5)
	rt= ResultsTable.getResultsTable()
	rt.reset() 
	analyzer = Analyzer(imp_,rt)
	imp_.setRoi(roi_)
	analyzer.measure()
	imp_.deleteRoi()
	imp_.getProcessor().subtract(rt.getValue ('Mean', 0))
	return

# Convert a freeline to a 2-pts PointRoi
def Line2PointRoi(roi_) :
	FPol = roi_.getInterpolatedPolygon(10, True) 
	points = PointRoi()
	points.addPoint(FPol.xpoints[0], FPol.ypoints[0])
	points.addPoint(FPol.xpoints[FPol.npoints-1], FPol.ypoints[FPol.npoints-1])
	points.setOptions("extra large cross label")
	return points

# Generic Dialog for threshold method
def thresholdImageUI(imp_):  
	IJ.setAutoThreshold(imp_, "Default dark")
	minThresholdValue = imp_.getProcessor().getMinThreshold()
	imp_.show()	
	gd = GenericDialog("Threshold")
	gd.addSlider("Min", 0, 65535, minThresholdValue)
	gd.addSlider("Max", 0, 65535, 65535)
	sliders = gd.getSliders()
	previewer = ThresholdPreviewer(imp_, sliders)  
	for slider in sliders :
		slider.addAdjustmentListener(previewer)  

	gd.showDialog()  

	if gd.wasOKed():  
		minThreshold = previewer.getMinThreshold()
		maxThreshold = previewer.getMaxThreshold()
		IJ.setRawThreshold(imp_, minThreshold, maxThreshold, None)
	else :
		IJ.setAutoThreshold(imp_, "MaxEntropy dark")
	imp_.hide()
	ImageConverter(imp_).convertToGray32()
	IJ.run(imp_, "NaN Background", "")	
	return


#Remove 0-value and saturated pixels and attribute the value "NAN" to them
def removeSaturatedPixels(imp_):
	imp2_= imp_.duplicate()
	IJ.setRawThreshold(imp2_, 1, 4094 ,None) # the camera is 12-bits dynamical range =[0,4095]
	IJ.run(imp2_, "Convert to Mask", "")
	IJ.run(imp2_, "Divide...", "value=255")
	ImageConverter(imp_).convertToGray32()
	IJ.run(imp_, "Gaussian Blur...", "sigma=2")
	imp_ = IC().run("Multiply create 32-bit", imp_, imp2_)
	IJ.setRawThreshold(imp_, 0, 4095 ,None)
	IJ.run(imp_, "NaN Background", "")	
	return

# Computation of the fret index
def CalculationFRETIndex(impD_,impA_, threshold):
	if threshold :
		impMask = impA_.duplicate()
		#Select the good threshold in the Acceptor image and threshold with NAN background
		#if not Autothreshold.MaxEntropy method is used
		thresholdImageUI(impMask)
		## Apply Mask to Donnor and Acceptor
		IC().run("Multiply", impD_,impMask)
		IC().run("Multiply", impA_,impMask)

	impD_ = IC().run("Add create 32-bit", impD_,impA_)#------- 1) image Donnor+Acceptor -> Denominator
	IJ.setRawThreshold(impD_,1,Float.MAX_VALUE, None) # remove 0-value pixels to exclude infinity value in the  divide calculation
	IJ.run(impD_, "NaN Background", "")	
	impA_ = IC().run("Divide create 32-bit", impA_, impD_)#------- 2) image Acceptor/(Donnor+Acceptor)
	IJ.setRawThreshold(impA_,0,1, None)
	IJ.run(impA_, "NaN Background", "")	
	IJ.run(impA_, "Multiply...", "value=100")
	IJ.run(impA_, "Enhance Contrast", "saturated=0.35")
	impA_.setDisplayRange(0, 100)
	IJ.run(impA_, "Fire", "")
	return impA_

			 	
# Label to ROIManager
def L2R(imp_, rm_):
	ip_ = imp_.duplicate().getProcessor()
	wand = Wand( ip_ )
	ip_.setColor(0)
	for y in range(imp_.height):
		for x in range(imp_.width):
			 if ( ip_.getPixel( x, y ) > 0.0 ):
			 	wand.autoOutline( x, y )
			 	if ( wand.npoints > 0 ) :
			 		roi = PolygonRoi(wand.xpoints, wand.ypoints, wand.npoints, Roi.TRACED_ROI)
					roi.setPosition( imp_.getCurrentSlice() )
					ip_.fill(roi)
					rm_.addRoi( roi )
	return

# check the value of the multipoint roi to know if it is necessary to merge the background with the function "mergeLabels"
def checkMerging(imp_,roi_):
	polyRoi = roi_.getFloatPolygon() 
	ip_= imp_.getProcessor()
	merging = False
	for i in range(polyRoi.npoints):
		merging |= ip_.getPixelValue(int(polyRoi.xpoints[i]), int(polyRoi.ypoints[i])) > .0
	return merging

#convert the 2 arrays of the polygon coordinates into Points2D (in order to use the vector product). Remove 2 points too close
def Polygon2Points(FPol) :
	p2p =[Point2D(FPol.xpoints[i], FPol.ypoints[i]) for i in range(FPol.npoints)]
	#remove all doublets
	return list(OrderedDict.fromkeys(p2p))

#---------------------------------------------------------------
#     -----------------       Start       -----------------
#---------------------------------------------------------------


#### STEP 1 : Preparation of data & analysis - Preprocessing

print "STEP 1 : Preparation of data & analysis - Preprocessing"
# clear the console automatically when not in headless mode
uiService.getDefaultUI().getConsolePane().clear()

#close Result Table if opened
if IJ.isResultsWindow() :
	IJ.run("Clear Results", "")
	tw = ResultsTable().getResultsWindow()
	tw.close()

#convert Files from #@ parameters to String 
imagefile = inputFile.getCanonicalPath()

# initialize the reader and get the OME metadata
reader = ImageReader()
omeMeta = MetadataTools.createOMEXMLMetadata()
reader.setMetadataStore(omeMeta)
reader.setId(imagefile)
sizeC = reader.getSizeC()
sizeT = reader.getSizeT()
sizeZ = reader.getSizeZ()
seriesCount = reader.getSeriesCount() #if multiple series 
reader.close()


if sizeC==1 or sizeT>1 or sizeZ>1 : #abort script if monochannel, multiframe or multislice
	IJ.error("nuclearFRET error", "The script requires a multichannel image with one slice and one frame ")
	sys.exit(0)

#select the serie if several
idxSerie=0
if seriesCount>1 :
	gui = GenericDialog("Select image serie")
	gui.addSlider("Image series: ", 1, seriesCount, 3)
	gui.showDialog()
	idxSerie = int(gui.getNextNumber()-1)

#Open a menu for selecting Donnor and Acceptor image indexes
print "Select Donnor, Acceptor images and the serie index to analyze"
idxDonnor, idxAcceptor, idProj= getImpIndexes(imagefile, sizeC, idxSerie)


#Extract Donnor and Acceptor images 
print "Process Donnor and Acceptor images" 
impDonnor = extractImpFromIndex(imagefile, idxDonnor, idxSerie)
impAcceptor = extractImpFromIndex(imagefile, idxAcceptor, idxSerie)
cal = impDonnor.getCalibration()
unit = cal.getUnit()
if unit =="micron" :
	unit ="um"


#### STEP 2 :  Subtract background and Calculate FRET index
print 'STEP 2 : FRET index image'

#Background subtraction
impProj = WindowManager.getImage(idProj)
removeSaturatedPixels(impDonnor)
removeSaturatedPixels(impAcceptor)

if backgroundSubtract:
	impDonnor.getProcessor().subtract(backgroundDonor)
	impAcceptor.getProcessor().subtract(backgroundAcceptor)
else :	
	print 'Select a ROI in the background'
	backROI = getBackgroundROI(impProj)
	subtractBackground(impDonnor,backROI)
	subtractBackground(impAcceptor,backROI)


#Create Folder for Donnor and Acceptor images & Save them
impFolder = createFolder(imagefile , idxSerie)
basename = os.path.basename(impFolder)

IJ.saveAs(impDonnor, "TIFF",os.path.join(impFolder, basename+"_c1.tif")) #save Donnor image
IJ.saveAs(impAcceptor, "TIFF",os.path.join(impFolder, basename+"_c2.tif")) #save Acceptor image


#close RoiManager --- All RM will not be visible
rm = RoiManager.getInstance()
if rm:
	rm.close()

impFRET = CalculationFRETIndex(impDonnor,impAcceptor, FRETtype == "Whole cell")
impFRET.setTitle(FRETTitle+".tif")
impFRET.setCalibration(cal)
IJ.saveAs(impFRET, "TIFF",os.path.join(impFolder, FRETTitle+".tif"))

if (FRETtype == "Whole cell") :
	impFRET.show()
	Analyzer.setMeasurements(Measurements.AREA+ Measurements.MEAN +Measurements.STD_DEV)
	Analyzer.setPrecision(5)
	rt= ResultsTable() 
	analyzer = Analyzer(impFRET,rt)
	analyzer.measure()
	rt.show("Mean FRET index (%)")
else : 
	#### STEP 3 :  Segmentation of nuclei  measurement
	print 'STEP 3 : Segmentation of nuclei'

	#Draw marker for MorphoLibJ Marker-controlled Watershed 
	print "Draw marker for Marker-controlled Watershed from MorphoLibJ library" 
	IJ.setTool("multipoint")
	impProj.show()
	IJ.run(impProj, "Enhance Contrast...", "saturated=0.3 equalize")
	rmMarker = RoiManager(showRoiManager)
	waitDialog = WaitForUserDialog("ROI Selection","Please select inside each object of interest using the point selection tool.")
	waitDialog.show()
	rmMarker.addRoi(impProj.getRoi())
	impProj.deleteRoi()
	if backgroundSubtract :
		backROI = getBackgroundROI(impProj)
	rmMarker.addRoi(Line2PointRoi(backROI))
	impProj.hide()
	

	bp = ByteProcessor(impAcceptor.width, impAcceptor.height)
	bp.setColor(Color.WHITE)
	for i in range(2):
		roi = rmMarker.getRoi(i)
		p = roi.getPolygon()
		for x,y in zip(p.xpoints, p.ypoints) :
			bp.fill(OvalRoi(x-2, y-2, 5, 5))
	impMarker = ImagePlus("Marker Image", bp) 
	impMarker = BinaryImages.componentsLabeling(impMarker, 8, 32)
	impLabel = Watershed.computeWatershed(impProj, impMarker, None, 8, True )
	lutName = CommonLabelMaps.JET.getLabel()
	lut = CommonLabelMaps.fromLabel(lutName).computeLut(255, True)    	
	impLabelRGB = LabelImages.labelToRgb(impLabel, lut ,Color.BLACK)
	impLabelRGB.show()

	#Merge regions of background and remove border label
	print "Merge regions of background and remove border label" 
	IJ.setTool("multipoint")
	waitDialog = WaitForUserDialog("ROI Selection","Please select background regions to merge using the point selection tool.")
	waitDialog.show()
	labelRoi = impLabelRGB.getRoi()
	rmMarker.addRoi(labelRoi)
	impLabelRGB.hide()
	if checkMerging(impLabel,labelRoi) :
		LabelImages.mergeLabels(impLabel, labelRoi, True)
	impLabelRGB.deleteRoi()
	LabelImages.removeBorderLabels(impLabel) 
	impLabelRGB = LabelImages.labelToRgb(impLabel, lut ,Color.WHITE)
	IJ.saveAs(impLabelRGB, "TIFF",os.path.join(impFolder, "LabelBordersRGB.tif"))
	rmMarker.runCommand("Deselect") # deselect ROIs to save them all
	rmMarker.runCommand("Save", os.path.join(impFolder, "RoiSet_Markers.zip")) #save the Rois


	#Select the periphery of the nuclei
	print "Select the periphery of the nuclei"
	rmNuclei = RoiManager(showRoiManager)
	rmNucleiOut = RoiManager(showRoiManager)
	L2R(impLabel, rmNuclei)
	count = rmNuclei.getCount() 
	rt = ResultsTable()

	rmContour = RoiManager(showRoiManager)
	ipFRET = impFRET.getProcessor()
	for idx in range(count): 
		roi0 = rmNuclei.getRoi(idx)
		FintPol = roi0.getInterpolatedPolygon(-1, True) 
		if FintPol.npoints > 10 : # exclude ROI with nb of coutour points < 10
			polyroi0 = PolygonRoi(FintPol, Roi.POLYGON)
			rmNuclei.addRoi(polyroi0)
			Fpts = Polygon2Points(FintPol)
			FptsSize = len(Fpts)
			for ipts in range(FptsSize) :
				Xpts = Fpts[ipts].getX()
				Ypts = Fpts[ipts].getY()	
				rt.incrementCounter()
				rt.addValue("IObject", idx)
				rt.addValue("IName", rmNuclei.getName(idx))
				rt.addValue("IContourPoints", ipts)
				rt.addValue("PointX", Xpts)
				rt.addValue("PointY", Ypts)
				rt.addValue(FRETTitle, ipFRET.getPixelValue(int(Xpts), int(Ypts)) )
		
			roiOUT = RoiEnlarger.enlarge(polyroi0, 2)
			rmNucleiOut.addRoi(roiOUT)
			roiIN = RoiEnlarger.enlarge(polyroi0, -1)
			notRoi = ShapeRoi(roiOUT).xor(ShapeRoi(roiIN))
			rmContour.addRoi(notRoi)
	rmNuclei.setSelectedIndexes(range(count))
	rmNuclei.runCommand(impLabel,"Delete")
	rmNuclei.runCommand("Deselect") # deselect ROIs to save them all
	rmNuclei.runCommand("Save", os.path.join(impFolder, "RoiSet_NucleiContour.zip")) #save the Contours
	
	for i in reversed(range(rt.getCounter())):
		if 	isnan(rt.getValue(FRETTitle, i)) or rt.getValue(FRETTitle, i) == 0:
			rt.deleteRow(i)


	IJ.run("Input/Output...", "jpeg=85 gif=-1 file=.csv save_column")
	rt.saveAs(os.path.join(impFolder,"ContourMeasurements.csv")) #save the measurement table

	#### STEP 4 :  FRET index of segmented nuclei
	print 'STEP 4 : FRET index of segmented nuclei'
	
	rt.reset()	
	Analyzer.setMeasurements (Measurements.AREA+ Measurements.MEAN +Measurements.STD_DEV + Measurements.SHAPE_DESCRIPTORS) 
	rmNucleiOut.runCommand(impFRET,"Measure")
	rt = Analyzer.getResultsTable()
	rt.saveAs(os.path.join(impFolder,"MeanFRETMeasurements.csv")) #save the measurement table
	IJ.renameResults("Mean FRET index (%)")
	rmContour.runCommand("Deselect") # deselect ROIs to save them all
	rmContour.runCommand("Save", os.path.join(impFolder, "RoiSet_NuclearBand.zip")) #save the nuclei band
	rmContour.setSelectedIndexes(range(rmContour.getCount()))
	rmContour.runCommand(impFRET,"Combine")
	allRoi = impFRET.getRoi()
	ipFRET.setValue(float('nan'))
	ipFRET.fillOutside(allRoi)	
	impFRET.deleteRoi()	
	impFRET.show()
	IJ.saveAs(impFRET, "TIFF",os.path.join(impFolder, "FRET index Nuclei.tif"))
	
print 'End'
