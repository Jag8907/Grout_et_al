
setBatchMode(true);

run("ImageJ2...", "scijavaio=true loglevel=INFO");

dir1= getArgument() ;

dir1= "D:/Samples/";

print("DIRECTORY: ");
print(dir1);


// markers can be replaced here 
cya = "[0-9a-zA-Z].*" + "hemato.*Colour_1" + ".*" ;
gra = "[0-9a-zA-Z].*" + "ker.*Colour_2" + ".*" ;
redi = "[0-9a-zA-Z].*" + "FAP" + ".*" ;
gre = "[0-9a-zA-Z].*" + "ADH1B" + ".*" ;
mag = "[0-9a-zA-Z].*" + "CD8" + ".*" ;
blu = "[0-9a-zA-Z].*" + "ASMA" + ".*" ;
yel = "[0-9a-zA-Z].*" + "CD34" + ".*" ;


dir2 = dir1 + "Registered/";
dir3 = dir1 +  "Color/";
dir4 = dir1 + "Overlay/";

dir2RGB = dir2 + "RGB/";
dir2hemato = dir2 + "hematoxylin/";
dir2HDAB = dir2 + "HDAB/";
dir2transf = dir2 + "transform/";
dir2align = dir2 + "aligned/";


// IMAGES ON WHICH TO APPLY PSEUDO COLORS 
run("ImageJ2...", "scijavaio=true loglevel=INFO");

//AVERAGING NUCLEI
//Images for averaging are ADH1B and CD34 - Could be change if necessary
print("AVERAGES");
list = getFileList(dir2align);
setBatchMode(true);
for (i=0; i<list.length; i++) {   
	showProgress(i+1, list.length);
	 if (matches(list[i], "[0-9a-zA-Z].*ADH1B.*Colour_1.*")) {
	 	toselect = dir2align +list[i] ;
	 	print(toselect) ;
	 	open(dir2align+list[i]);
	 	open(dir2align+list[i]);}
	 else if (matches(list[i], "[0-9a-zA-Z].*CD34.*Colour_1.*")) {
	 	toselect = dir2align+list[i] ;
	 	print(toselect) ;
	 	open(dir2align+list[i]);}
	  	
}
run("Concatenate...", "all_open title=hemato_Colour_1_2");
run("Z Project...", "projection=[Average Intensity]");
close("hemato_Colour_1_2");

list = getFileList(dir2align);
setBatchMode(true);
for (i=0; i<list.length; i++) {   
	showProgress(i+1, list.length);
	if (endsWith(list[i], "Colour_2).tif")) {
	 	open(dir2align+list[i]); }
}


ids=newArray(nImages); 
stringmerge = "";


setOption("ExpandableArrays", true);
vectColor = newArray;

print("INITIALIZATION OK");
	
for (i=0;i<nImages;i++) {
	selectImage(i+1); 
	title = getTitle;
	print(title);
	ids[i]=getImageID;
	//Line to use if we want to save  not tresholded image in the colour folder (idem for each marker)
	//saveAs("tiff", dir3 + title );
	run("Invert")  ;
	if (i+1==1){
		imwidth = getWidth() ;
		imheigth = getHeight();
		overlap = 0.025 * imwidth ;
		print("dimensions");
		print(imwidth);
		print(imheigth);
		}
	makeRectangle(overlap, overlap, imwidth - 2* overlap, imheigth - 2* overlap);
	run("Crop");
	changeValues(0,0, NaN);
	changeValues(255,255, NaN);
	print("newdim");
	print(getWidth);
	print(getHeight);
	k = i+1;
	if (startsWith(title, "no_stain") == false) {stringmerge = stringmerge + "image"+k+"=[" + title + "] " ;}
	print(title);
	if (matches(title, redi)==1){
		run("8-bit");
		run("Red"); 
		//saveAs("tiff", dir3 + title + "_Red" );
		setAutoThreshold("MaxEntropy dark");
		getThreshold(lower, upper);
		lower2 = lower - 10 ;
		setMinAndMax(lower2, 200);
		run("Apply LUT");
		imageTitle1 = getTitle();
		if (startsWith(title, "no_stain") == false) {vectColor[i] = "Red";}
		}
	else if (matches(title, gre)==1) {
		run("8-bit");
		run("Green");
		//saveAs("tiff", dir3 + title  + " _Green"); 
		setAutoThreshold("MaxEntropy dark");
		getThreshold(lower, upper);
		lower2 = lower - 5 ;
		setMinAndMax(lower2, 200);
		run("Apply LUT");
		imageTitle2 = getTitle();
		if (startsWith(title, "no_stain") == false) {vectColor[i] = "Green";}
		}
	else if (matches(title, blu)==1) {
		run("8-bit");
		run("Blue");
		//saveAs("tiff", dir3 + title  + " _Blue");
		setAutoThreshold("MaxEntropy dark");
		getThreshold(lower, upper);
		lower2 = lower - 5 ;
		setMinAndMax(lower2, 200);
		run("Apply LUT");
		imageTitle3 = getTitle();
		if (startsWith(title, "no_stain") == false) {vectColor[i] = "Blue";}
		}
	else if (matches(title, gra)==1) {
		run("8-bit");
		run("Grays");
		//saveAs("tiff", dir3 + title + "_Grays");
		setAutoThreshold("Otsu dark");
		getThreshold(lower, upper);
		setMinAndMax(lower, 200);
		run("Apply LUT");
		imageTitle4 = getTitle(); 
		if (startsWith(title, "no_stain") == false) {vectColor[i] = "Grays";}  
		}
	else if (matches(title,  cya)==1) {
		run("8-bit");
		run("Cyan");
		//saveAs("tiff", dir3 + title + "_Cyan");
		setAutoThreshold("MaxEntropy dark");
		getThreshold(lower, upper);
		setMinAndMax(lower, 200);
		run("Apply LUT");
		imageTitle5 = getTitle(); 
		print(imageTitle5);
		if (startsWith(title, "no_stain") == false) {vectColor[i] = "Cyan";} 
		} 
	else if (matches(title, mag)==1) {
		run("8-bit");
		run("Magenta");
		//saveAs("tiff", dir3 + title + "_Magenta" );
		setAutoThreshold("MaxEntropy dark");
		getThreshold(lower, upper);
		setMinAndMax(lower, 200);
		run("Apply LUT");
		imageTitle6 = getTitle();
		if (startsWith(title, "no_stain") == false) {vectColor[i] = "Magenta";} 
		}
	else if (matches(title, yel)==1) {
		run("8-bit");
		run("Yellow");
		//saveAs("tiff", dir3 + title + "_Magenta" );
		setAutoThreshold("Otsu dark");
		getThreshold(lower, upper);
		lower2 = lower - 10 ;
		setMinAndMax(lower2, 200);
		run("Apply LUT");
		imageTitle6 = getTitle();
		if (startsWith(title, "no_stain") == false) {vectColor[i] = "Yellow";} 
		}
	else {
		run("8-bit");
		run("Yellow");
		//saveAs("tiff", dir3 + title + "_Yellow" ); 
		setAutoThreshold("MaxEntropy dark");
		getThreshold(lower, upper); 
		lower2 = lower - 5 ;
		setMinAndMax(lower2, 200);
		run("Apply LUT");
		imageTitle7 = getTitle();
		if (startsWith(title, "no_stain") == false) {vectColor[i] = "Yellow";} 
		}

		//saveAs("tiff", dir3 + title + "_modif" );
}

print(stringmerge);
print(stringmerge + "image"+nImages+"=[-- None --]");


nchannels = nImages ;
nchannelsplus = nchannels + 1;
print(nchannelsplus);

run("Concatenate...", "  title=Stack "+ stringmerge + "image"+nchannelsplus+"=[-- None --]");

run("Stack to Hyperstack...", "order=xyczt(default) channels="+nchannels+" slices=1 frames=1 display=Composite");

for(c=0;c<lengthOf(vectColor);c++) {
    Stack.setChannel(c+1);
    currentcol = vectColor[c];
    
    if (startsWith(toString(vectColor[c]), "Red") == true) {run("Red");}
    else if (matches(toString(vectColor[c]), "Blue") == true) {run("Blue");}
    else if (matches(toString(vectColor[c]), "Green") == true) {
    	run("Green");
    	
    	}
    else if (matches(toString(vectColor[c]), "Grays") == true) {run("Grays");}
    else if (matches(toString(vectColor[c]), "Magenta") == true) {run("Magenta");}
    else if (matches(toString(vectColor[c]), "Cyan") == true) {run("Cyan");}
    else {run("Yellow");}
    
    print(currentcol);
    
}

print("saving");

ending = indexOf(dir1, "ROIsTif") ;
//If the name of the folder is not 3 letters -> for ex : 62P4 put -5 ; eg 24_ROIsTiF => c24_composite.ome.tif     c is for compressed 
beginning = ending -5 ;
title = substring(dir1, beginning, ending);
title = title + "composite";
saveAs("tiff", dir4 + title );
run("Bio-Formats Exporter", "save="+ dir4  + "c" + title  + ".ome.tif compression=LZW") ;

close("*") ;
run("Collect Garbage");
call("java.lang.System.gc");

print("END") ; 
File.delete(dir2RGB);
File.delete(dir3);
run("ImageJ2...", "scijavaio=false loglevel=INFO");

eval("script", "System.exit(0);"); 

