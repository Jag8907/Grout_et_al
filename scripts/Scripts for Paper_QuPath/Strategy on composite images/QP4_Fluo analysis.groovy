
/*
//To perform manually: draw ROI, set to other, change channel names if needed. Duplicate. Definition of all the classifiers (TumorvsStroma ; ASMAclassif ; FAPvsADH1B).
setPixelSizeMicrons(0.4414, 0.4414)
setChannelNames('nuclei', 'keratin', 'ADH1B', 'CD34', 'CD8', 'FAP', 'MYH11', 'PDPN')
*/

def PathExport = "D:/Samples/results/"

// Get areas per tiles ====================================================================================================
selectObjects { p -> p.getPathClass() == getPathClass("Other") && p.isAnnotation() }
runPlugin('qupath.lib.algorithms.TilerPlugin', '{"tileSizeMicrons": 500.0,  "trimToROI": true,  "makeAnnotations": true,  "removeParentAnnotation": true}');
def annotationstiles = getAnnotationObjects()
annotationstiles.each {it.setPathClass(getPathClass('Other'))}


selectObjects { p -> p.getPathClass() != getPathClass("None") && p.isAnnotation() }
createAnnotationsFromPixelClassifier("TumorvsStroma", 100.0, 1000.0, "DELETE_EXISTING", "SELECT_NEW")
createAnnotationsFromPixelClassifier("ASMAclassif", 5.0, 30.0)

//================================================================================================================================================================
// Get the list of all images in the current project
def project = getProject()
def imagesToExport = getCurrentImageData();
//def imagesToExport = project.getImageList()
def CurrentName = GeneralTools.getNameWithoutExtension(imagesToExport.getServer().getMetadata().getName())
print("Current image name: " + CurrentName)



// Choose your *full* output path
def outputPath1 = PathExport + CurrentName + "_annotation_results_tiles.txt"
def outputFile = new File(outputPath1)
saveAnnotationMeasurements(outputPath1)
print "annotation saved!"



// Global analysis===========================================================================================================


selectObjectsByClassification("Tumor");
mergeSelectedAnnotations()

selectObjectsByClassification("Stroma");
mergeSelectedAnnotations()

selectObjectsByClassification("ASMA");
clearSelectedObjects();


//selectObjectsByClassification("Necrosis");
//mergeSelectedAnnotations()

//selectObjectsByClassification("Alveoli");
//mergeSelectedAnnotations()


selectObjectsByClassification("Other");
clearSelectedObjects();

print "annotation merged"

selectAnnotations();
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons": 2.0,  "region": "ROI",  "tileSizeMicrons": 25.0,  "channel1": false,  "channel2": false,  "channel3": false,  "channel4": true,  "channel5": false,  "channel6": false,  "channel7": false,  "doMean": true,  "doStdDev": false,  "doMinMax": false,  "doMedian": false,  "doHaralick": false,  "haralickMin": NaN,  "haralickMax": NaN,  "haralickDistance": 1,  "haralickBins": 32}');


selectObjects { p -> p.getPathClass() != getPathClass("None") && p.isAnnotation() }
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage": "nuclei",  "requestedPixelSizeMicrons": 0.4,  "backgroundRadiusMicrons": 8.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 10.0,  "maxAreaMicrons": 400.0,  "threshold": 15.0,  "watershedPostProcess": true,  "cellExpansionMicrons": 3.0,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

print "Cell detection done"


detectionToAnnotationDistances(true)

print "distances done"

//Extend annotation=========================================================================================================



import qupath.lib.gui.tools.MeasurementExporter
import qupath.lib.objects.PathCellObject


def myClass1 = getPathClass('20µm')
def myClass2 = getPathClass('10µm')



import qupath.lib.gui.measure.ObservableMeasurementTableData
import qupath.lib.objects.PathCellObject

// Get observable measurement table for the current image
def annotations = getAnnotationObjects()
def ob = new ObservableMeasurementTableData()
ob.setImageData(getCurrentImageData(),  annotations)


mini = 0

annotations.each {
    area= ob.getNumericValue(it, "Area µm^2")
    if (it.getPathClass() == getPathClass('Tumor')){mini = area}
   
}

print(mini)


//20 µm
selectObjectsByClassification("Tumor");
runPlugin('qupath.lib.plugins.objects.DilateAnnotationPlugin', '{"radiusMicrons": 20.0,  "lineCap": "Round",  "removeInterior": false,  "constrainToParent": true}');

def annotations1 = getAnnotationObjects()
def ob1 = new ObservableMeasurementTableData()
ob.setImageData(getCurrentImageData(),  annotations1)

def matrixArea1 = []
annotations1.each { 
    area= ob.getNumericValue(it, "Area µm^2")
    print(area)
    matrixArea1 << area
    if (it.getPathClass() == getPathClass('Tumor') && area > mini) {it.setPathClass(getPathClass('20µm'))}  
}



//10 µm
selectObjectsByClassification("Tumor");
runPlugin('qupath.lib.plugins.objects.DilateAnnotationPlugin', '{"radiusMicrons": 10.0,  "lineCap": "Round",  "removeInterior": false,  "constrainToParent": true}');

def annotations2 = getAnnotationObjects()
def ob2 = new ObservableMeasurementTableData()
ob.setImageData(getCurrentImageData(),  annotations2)

def matrixArea2 = []
annotations2.each { 
    area= ob.getNumericValue(it, "Area µm^2")
    print(area)
    matrixArea2 << area
    if (it.getPathClass() == getPathClass('Tumor') && area > mini) {it.setPathClass(getPathClass('8µm'))}  
}

// ADH1B & FAP =======================================================

selectObjectsByClassification("Stroma");
createAnnotationsFromPixelClassifier("FAPvsADH1B", 10.0, 30.0)
createAnnotationsFromPixelClassifier("ASMAclassif", 5.0, 30.0)
selectAnnotations();
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons": 2.0,  "region": "ROI",  "tileSizeMicrons": 25.0,  "channel1": false,  "channel2": false,  "channel3": false,  "channel4": true,  "channel5": false,  "channel6": false,  "channel7": false,  "doMean": true,  "doStdDev": false,  "doMinMax": false,  "doMedian": false,  "doHaralick": false,  "haralickMin": NaN,  "haralickMax": NaN,  "haralickDistance": 1,  "haralickBins": 32}');


// =========================================================================
//SAVING distance to tumor border

def outputPath2 = PathExport + CurrentName + "_annotation_results_total.txt"
//saveDetectionMeasurements(outputPath2 )
saveAnnotationMeasurements(outputPath2)

def outputPath3 = PathExport + CurrentName + "_detection_distances.txt"
//saveDetectionMeasurements(outputPath3, "Parent", "Centroid X µm", "Centroid Y µm" )
saveDetectionMeasurements(outputPath3 )

Thread.sleep(100)

// Try to reclaim whatever memory we can, including emptying the tile cache

javafx.application.Platform.runLater {

    getCurrentViewer().getImageRegionStore().cache.clear()

    System.gc()

}

Thread.sleep(100)
