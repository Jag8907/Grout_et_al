setBatchMode(true);

run("ImageJ2...", "scijavaio=true loglevel=INFO");

dir1= getArgument() ;

dir1= "D:/Samples/";

print("DIRECTORY: ");
print(dir1);

list = getFileList(dir1);


dir2 = dir1 + "Registered/";
dir3 = dir1 +  "Color/";
dir4 = dir1 + "Overlay/";

dir2RGB = dir2 + "RGB/";
dir2hemato = dir2 + "hematoxylin/";
dir2HDAB = dir2 + "HDAB/";
dir2transf = dir2 + "transform/";
dir2align = dir2 + "aligned/";


setBatchMode(true);

runMacro("D:/alignment_hemato_local.py",  dir2);

close("*");
run("Collect Garbage");
call("java.lang.System.gc");

eval("script", "System.exit(0);"); 

