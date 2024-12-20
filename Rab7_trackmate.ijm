//E. Gutierrez-Galindo et al, University of Stuttgart (2024)

//This macro was used to find and segment Rab7-vesicles using maximum intensity projections (MIP), and save them as ROIs. It uses the plugin "Trackmate" (Ershov et al., 2022).
//IMPORTANT: These saved ROIs (Rab7 vesicles) are then measured (area, count, intensity...) using the macro "Rab7_trackmate_measureRIOs.ijm"
//Requires: 1) MIPs ("input_path"); 2) Cell ROIs (cell outlines) saved in "roi_path" - this can be saved manually or using Save_ROIs macro


/////////////////
///   MAIN   ///
////////////////

//Initiate Variables
input_path = getDirectory("folder with input images");
roi_path = getDirectory("folder where ROIs are saved");
file_list = getFileList(input_path);
output_path = getDirectory("folder where analysed images are saved")
ROI_path = getDirectory("folder where ROIs for TM vesicles are saved")


for(P = 0; P < file_list.length; P++) {print(file_list[P]);
	if (endsWith(file_list[P], ".czi")) {
	
		run("Bio-Formats Importer", "open=[" + input_path + file_list[P] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension;
		title_raw = getTitle();
		title_wo_MIP = substring(title, 0, lengthOf(title)-4);
		load_ROIs(title_wo_MIP, roi_path);
		cell_no = make_Cell_pics(title_raw);
		for (i = 0; i < cell_no; i++) {
			trackmate(title + "_cell" + i);
		}
		
		run("Close All");
		roiManager("reset");

	}
}


//////////////////////
///   FUNCTIONS   ///
/////////////////////

function make_Cell_pics(sel_wind) { 
// function to make duplicates for each cell
	cells = roiManager("Count");
	selectWindow(sel_wind);
	for (j = 0; j < cells; j++){
		selectWindow(sel_wind);
		roiManager("Select", j);
		run("Duplicate...", "title=[" + title + "_cell"+ j + "] duplicate channels=2"); //select Rab7 channel
		run("Clear Outside", "stack");

	}
	roiManager("reset");
	return cells;

}

function trackmate(image) {

	run("Set Measurements...", "area mean integrated display redirect=None decimal=3");
	selectWindow(image);
	
	run("Duplicate...", "title=[" + image + "_sub-bgr]");
	run("Subtract Background...", "rolling=1 sliding");
	run("Duplicate...", "title=[" + image + "_TM]");
	
	
	//This step is not automated (needs to be done manually for each cell, and then press "ok" to continue to the next cell)
	run("TrackMate"); //1.8, 3 (settings used for all cells)
	waitForUser("Press OK after TM");
	run("Flatten");
	
	
	saveAs("Jpeg", output_path + image + "_vesicles" + ".jpeg");
	
	
	selectWindow(image + "_sub-bgr");
	ROI_amount = roiManager("count");
	
	if (ROI_amount==0){
			continue;
		}
		else {roiManager("Save", ROI_path + image + ".zip");
		}
	

	roiManager("reset");
	
}	
	

function load_ROIs(title, path_roi) {
	roiManager("Open", path_roi + title + ".zip");
}

		
