double downsample = 1

def server = getCurrentServer()

def baseImageName = getProjectEntry().getImageName()



path = buildFilePath(PROJECT_BASE_DIR)
print 'path'
print path
int last = path.trim().length()
currentSample = path.substring(last -4, last);

path2 = buildFilePath(PROJECT_BASE_DIR, currentSample + '_ROIsTif_DS' + downsample)
mkdirs(path2)

def pathROI = "D:/Samples/"
 
 
getAnnotationObjects().eachWithIndex{it, x->
        println("Working on: "+it)
    def roi = it.getROI()
    def requestROI = RegionRequest.createInstance(server.getPath(), downsample, roi)
    def imagesToExport = getCurrentImageData();
    def CurrentName = GeneralTools.getNameWithoutExtension(imagesToExport.getServer().getMetadata().getName())
    currentImagePath = buildFilePath(PROJECT_BASE_DIR, currentSample + '_ROIsTif_DS' + downsample ,CurrentName +'_'+ x + '.ome.tif')
    writeImageRegion(server, requestROI, currentImagePath)
    
    currentImagePath2 = buildFilePath(pathROI, currentSample + '_ROIsTif_DS' + downsample ,CurrentName +'_'+ x + '.ome.tif')
    writeImageRegion(server, requestROI, currentImagePath2)
}
print 'DONE!'