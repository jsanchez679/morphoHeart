/*
 *  pA_ECandDx3.ijm
 *  This macro will:
 *  - Enhance contrast using 0.3% saturated pixels and normalizing the histogram 
 *  - Enhance contrast equalizing the histogram 
 *  - Despeckle (x3) 
 *  - Crop both channels using the MIP as a guidance 
 *  - Save each cropped channel individually _EDC.tif
 *  
 *  Created by: 
 *  Juliana Sanchez-Posada v1.0 March 30/20
 *  Corrected bugs:
 *  	
*/

macro "pA_EC_Dx3_Crop" {
	print("\\Clear")
	print("NEW RUN");
	
	// Get images source directory 
	dirSource = getDirectory("Choose the source directory where the individual channels are saved... ");
	print("dirSource: "+dirSource);

	Dialog.create("Get filename (LSXX_FXX_V_XX_XXXX_EDC):");
	Dialog.addString("LS (XX):", "");
	Dialog.addString("F (XX):", "");
	Dialog.addString("V/D (X):", "V");
	Dialog.addString("SS/DS (XX):", "DS");
	Dialog.addString("Time (XXXX):", "");
	
	Dialog.show();
	numLS = Dialog.getString();
	numF = Dialog.getString();
	VoD = Dialog.getString();
	SSoDS = Dialog.getString();
	timeT = Dialog.getString();
	
	filename = "LS"+numLS+"_F"+numF+"_"+VoD+"_"+SSoDS+"_"+timeT;
	print("filename: "+filename);

	//Channel 1 - After Arivis - Change name from ch1 to ch0
	ch1_filename = dirSource + filename+"_ch0";
	print("ch1_filename: "+ch1_filename);
	
	//Channel 2 - After Arivis - Change name from ch2 to ch1
	ch2_filename = dirSource + filename+"_ch1";
	print("ch2_filename: "+ch2_filename);

	run("ROI Manager...");
		
	//Open channel 1
	waitForUser("Open channel 1 and click OK when ready.");
	ch1_ID = getImageID();
	selectImage(ch1_ID);
	fECandDx3();
	ch1ECD_ID = getImageID();
	ch1ECD_tt = getTitle();

	//Open channel 2
	waitForUser("Open channel 2 and click OK when ready.");
	ch2_ID = getImageID(); 
	fECandDx3();
	ch2ECD_ID = getImageID();
	ch2ECD_tt = getTitle();

	//Merge channels
	Label_Comp = "c2="+ch1ECD_tt+" c6="+ch2ECD_tt+" create keep";
	run("Merge Channels...", Label_Comp);
	//Max Intensity Projection
	run("Z Project...", "projection=[Max Intensity]");
	MAX_ID = getImageID();

	// Set line width to 10
	run("Line Width...", "line=10");
	run("Colors...", "foreground=black background=black selection=yellow");

	//Crop image
	//Select MAX image to draw square to crop
	selectImage(MAX_ID);
	setTool("rectangle");
	waitForUser("Draw square -Shift- that encompasses the whole heart and click OK when ready");
	roiManager("Add");

	selectImage(ch1_ID);
	roiManager("Select", 0);
	run("Crop");
	setTool("rectangle");
	waitForUser("Draw rectangle that encloses the whole image and avoids open contours and click OK when ready");
	if (selectionType() == 0) {
		roiManager("Add");
		run("Draw", "stack"); // same as Ctrl+D
		print("Black rectangle drawn in ch1");
	}
	makeRectangle(0, 0, 1, 1);
	saveAs("tif", ch1_filename+"_EDC");
	ch1ECD_ID = getImageID();
	ch1ECD_tt = getTitle();
	
	selectImage(ch2_ID);
	roiManager("Select", 0);
	run("Crop");
	roiManager("Select", 1);
	run("Draw", "stack"); // same as Ctrl+D
	print("Black rectangle drawn in ch2");
	makeRectangle(0, 0, 1, 1);
	saveAs("tif", ch2_filename+"_EDC");
	ch2ECD_ID = getImageID();
	ch2ECD_tt = getTitle();

	waitForUser("Check stacks and click OK when ready to continue");
	selectImage(ch1ECD_ID);
	masking();
	ch1mask_ID = getImageID();
	ch1mask_tt = getTitle();
	saveAs("tif", ch1_filename+"_mask");

	selectImage(ch2ECD_ID);
	masking();
	ch2mask_ID = getImageID();
	ch2mask_tt = getTitle();
	saveAs("tif", ch2_filename+"_mask");

	waitForUser("Check stacks and click OK when ready to close");
	close("*");
	roiManager("reset")
	
	close("*");
	close("\\Others");
	
	print("Images have been closed");
	print("DONE: pA_EC_Dx3_Crop_Mask");
}

function fECandDx3() {
	run("Enhance Contrast...", "saturated=0.3 normalize process_all");
	run("Despeckle", "stack");
	run("Despeckle", "stack");
	run("Despeckle", "stack");
}

function masking() {
	run("Morphological Filters (3D)", "operation=Laplacian element=Cube x-radius=2 y-radius=2 z-radius=2");
	run("Invert", "stack");
	run("Threshold...");
	setAutoThreshold("Li dark");
	waitForUser("Set threshold (DO NOT CLICK APPLY) and click OK when ready");
	run("Convert to Mask", "method=Li background=Dark black");
	run("Invert", "stack");
	run("Morphological Filters (3D)", "operation=Dilation element=Cube x-radius=1 y-radius=1 z-radius=1");
}

