import os
from ij import IJ, ImagePlus
from ij.gui import GenericDialog
from register_virtual_stack import Register_Virtual_Stack_MT
#import register_virtual_stack.Transform_Virtual_Stack_MT
from register_virtual_stack import Transform_Virtual_Stack_MT
import fnmatch
from xml.etree import cElementTree as ET
import mpicbg.trakem2.transform.CoordinateTransform
import mpicbg.trakem2.transform.CoordinateTransformList
import gc

dire = getArgument()
print(dire)

# source directory
source_dir = dire+"hematoxylin/"
source_dir2 = dire+"HDAB/"
# output directory
target_dir = dire + "aligned/"
# transforms directory
transf_dir = dire+"transform/"
transf_dir2 = dire+"transform2/"

use_shrinking_constraint = 0
 

p = Register_Virtual_Stack_MT.Param()
# SIFT parameters
p.sift.maxOctaveSize = 4500 # ref 1024 for smaller images
#p.sift.fdSize = 5 # ref 4
#p.sift.initialSigma = 1.5 #ref 1.6
#p.sift.steps = 5    #ref 3
#p.maxEpsilon = 15
p.featuresModelIndex = 3 


# The "inlier ratio":
p.minInlierRatio = 0.1
 

source_dir = dire+"hematoxylin/"
for fileHemato in os.listdir(source_dir):
	if fnmatch.fnmatch(fileHemato, '*ker*'):
		reference_name = fileHemato
		print(reference_name)

print( "alignement hemato")
Register_Virtual_Stack_MT.exec(source_dir, target_dir, transf_dir, reference_name, p, use_shrinking_constraint)

IJ.run("Close All", "");

print("step python DONE") 
