setBatchMode(true);

run("ImageJ2...", "scijavaio=true loglevel=INFO");

dir1= getArgument() ;

dir1="D:/samples/";

print("DIRECTORY: ");
print(dir1);


list = getFileList(dir1);


dir2 = dir1 + "Registered/";
dir3 = dir1 +  "Color/";
dir4 = dir1 + "Overlay/";
File.makeDirectory(dir2);
File.makeDirectory(dir3);
File.makeDirectory(dir4);
dir2RGB = dir2 + "RGB/";
dir2hemato = dir2 + "hematoxylin/";
dir2HDAB = dir2 + "HDAB/";
dir2transf = dir2 + "transform/";
dir2align = dir2 + "aligned/";
File.makeDirectory(dir2RGB);
File.makeDirectory(dir2hemato);
File.makeDirectory(dir2HDAB);
File.makeDirectory(dir2transf);
File.makeDirectory(dir2align);



for (i=0; i<list.length; i++) {
     showProgress(i+1, list.length);
     print(list[i]);
      if (endsWith(list[i], ".tif")) {
          toselect = dir1+list[i] ;
          print(toselect) ; }
     run("Bio-Formats Importer", "open=[" + toselect + "] autoscalecolor_mode=Default view=Hyperstack stack_order=XYCZT");
     title = getTitle;
     run("RGB Color");
     selectWindow(title);
     close();
     run("Collect Garbage");
     call("java.lang.System.gc");
     ids=newArray(nImages);
     title = getTitle;
     print(title);
     // Remove black points - Change it/remove it if necessary/not needed
     changeValues(0x000000,0x102314,0xfefefe);
     print("new deconvo");
     run("Colour Deconvolution", "vectors=[H DAB]");
     print("SAVING hemato hdab");
     close("*Colour_3*");
     close("*(RGB).tif");
     title = getTitle;
     print(title);
     if (matches(title, "[0-9a-zA-Z].*Colour_2.*")==1){saveAs("tiff",dir2HDAB + title ) ;}
     close("*Colour_2*");
     title = getTitle;
     print(title);
     if (matches(title, "[0-9a-zA-Z].*Colour_1.*")==1){saveAs("tiff",dir2hemato + title ) ; }
     close("*");
     run("Collect Garbage");
     call("java.lang.System.gc");

     }


close("*");
run("Collect Garbage");
call("java.lang.System.gc");
eval("script", "System.exit(0);"); 
