run("Set Measurements...", "centroid redirect=None decimal=3");
setBatchMode(true);
run("To ROI Manager");
nuc_num=roiManager("count");

WinName=getTitle();
idx=indexOf(WinName ".svs");
SlideName=substring(WinName,0, idx)+"_J3";


save_path="//kuehlapis/scratch/shiva/jokhun/Project 4 - Tissue Microarray/H&E from BioMax/IndiNuc(New)/";

for (i=0; i<nuc_num; i++) {
	 
	roiManager("Select", i);

	roiManager("Measure");
	X=getResult("X",0);
	Y=getResult("Y",0);
	
	//if(d2s(X,1) !="NaN"){

	file_name=SlideName+"_"+(i+1)+"_"+X+"_"+Y;
	run("Clear Results");

	run("Duplicate...", " ");
	setBackgroundColor(0, 0, 0);
	run("Clear Outside");

	saveAs("Tiff", save_path+file_name);
	close();
	//}
	//else{
	//run("Clear Results");
	//}
	

}


selectWindow(WinName);
close();
