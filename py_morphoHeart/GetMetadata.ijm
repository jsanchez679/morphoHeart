/*
 *  GetMetadata.ijm
 *  This macro will get all the Metadata from the .czi, .nd2 files inside a folder and save it as a .txt file within that folder. 
 *  
 *  Juliana Sanchez-Posada Version Feb 12, 2020
 *   Corrected bugs:
 *  	XX-XXX-XX 	
 *  	
*/

macro "GetMetadata" {
	print("\\Clear")
	print("NEW RUN");
	run("Bio-Formats Macro Extensions");
	//Dialog to select the extension type of the files
	Dialog.create("Select the extension of the files you want to get the Metadata from");
	// If you want to add other type of extension, just add it to the list. NOTE: The extension has to be supported by Bioformats PlugIn in order to work.
	Dialog.addChoice("File extension:", newArray(".czi", ".nd2", ".TIF",".tif")); 
	Dialog.show();

	fileType = Dialog.getChoice();
	print("fileType: " + fileType);
	
	//Select the source directory
	setOption("JFileChooser", true); 
	dirSource = getDirectory("Go to the directory where the images you want to get the Metadata from are saved... ");
	run("Bio-Formats Macro Extensions");

	//Get the file list on the source directory
	list = getFileList(dirSource);
	print("n:"+list.length);
	numIm2An = 0; 
	
	//Creating txt file
	DataCode = getString("Enter the code for the data to get metadata from (e.g. LSXX, ASXX):","")
	MetaData_FileName= "METADATA_RAW_"+DataCode+".txt";
	MetaData_File=File.open(dirSource+MetaData_FileName);
	print(MetaData_File,dirSource+"\t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t \t ")
	print(MetaData_File,"ImageTitle\tFileName\tCreationDate\tTimeSeries\tLaser1\tLaserPower1\tLaser2\tLaserPower2\tFilter1\tFilter2\tFilter3\tFilter4\tCamFrameHeight\tCamFrameWidth\tScalingX\tScalingY\tScalingZ\tDimensionX\tDimensionY\tDimensionZ");

	for (i=0; i<list.length; i++) {
		filename = dirSource + list[i];
		print("filename: "+filename);
		
			if (endsWith(filename,fileType)) {
				numIm2An = numIm2An +1;
				//Set batch mode
				setBatchMode(true);
				print("Batch mode true");
				//Open file using Bioformats if it is a Z-Stack or normal open macro if it is a single image
				txt2run = "open=["+filename+"] autoscale color_mode=Composite display_metadata rois_import=[ROI manager] view=[Metadata only]  stack_order=Default";
				//run("Bio-Formats Importer", "open=[dir] autoscale color_mode=Composite display_metadata rois_import=[ROI manager] view=[Metadata only] stack_order=Default");
				run("Bio-Formats Importer", txt2run);//"open="+filename+" autoscale color_mode=Composite display_metadata rois_import=[ROI manager] view=[Metadata only]  stack_order=Default");	
				//Get Dimensions
				//Stack.getDimensions(widthNum, heightNum, channelsNum, slicesNum, framesNum);
				imageTitle = list[i];
				print("Image Title:" + imageTitle);
				showStatus("Getting Metadata from File: " + imageTitle);
				print("Getting Metadata from File: " + imageTitle);
				
				Ext.setId(filename);
				Ext.getMetadataValue("Information|Document|Name #1", md_FileName);
				
				Ext.getMetadataValue("Information|Document|CreationDate #1", md_CreationDate);
				
				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|TimeSeries #1", md_TimeSeries);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|Attenuator|Laser #1", md_Laser1);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|Laser|LaserPower #1", md_LaserPower1);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|Attenuator|Laser #2", md_Laser2);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|Laser|LaserPower #2", md_LaserPower2);

				Ext.getMetadataValue("Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|BeamSplitter|Filter #1", md_Filter1);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|BeamSplitter|Filter #2", md_Filter2);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|BeamSplitter|Filter #3", md_Filter3);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|BeamSplitter|Filter #4", md_Filter4);

				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|CameraFrameHeight #1", md_CamFrameHeight);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|CameraFrameWidth #1", md_CamFrameWidth);
				
				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingX #1", md_ScalingX);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingY #1", md_ScalingY);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|ScalingZ #1", md_ScalingZ);


				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|DimensionX #1", md_DimensionX);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|DimensionY #1", md_DimensionY);
				Ext.getMetadataValue("Experiment|AcquisitionBlock|AcquisitionModeSetup|DimensionZ #1", md_DimensionZ);

				print(MetaData_File,imageTitle+"\t"+md_FileName+"\t"+md_CreationDate+"\t"+md_TimeSeries+"\t"+md_Laser1+"\t"+md_LaserPower1+"\t"+md_Laser2+"\t"+md_LaserPower2+"\t"+md_Filter1+"\t"+md_Filter2+"\t"+md_Filter3+"\t"+md_Filter4+"\t"+md_CamFrameHeight+"\t"+md_CamFrameWidth+"\t"+md_ScalingX+"\t"+md_ScalingY+"\t"+md_ScalingZ+"\t"+md_DimensionX+"\t"+md_DimensionY+"\t"+md_DimensionZ);

				Ext.close()
				close("*");
			}}

	listFinal = getList("window.titles");
	for (i=0; i<listFinal.length; i++){
     	winame = listFinal[i];
      	selectWindow(winame);
     	run("Close");
     }
      
//run("Close");
print("FINISHED");
showStatus("finished");
}
     
