//E. Gutierrez-Galindo et al, University of Stuttgart (2024)

//This macro was used to measure area, count, intensity of Rab7-vesicles segmented with Trackmate
//IMPORTANT: First use "Rab7_trackmate.ijm" to segment Rab7 vesicles
//Requires: 1) MIPs ("input_path"); 2) Cell ROIs (cell outlines) saved in "roi_path" - this can be saved manually or using Save_ROIs macro; 3) Trackmate ROIs saved (ROI_path2).


/////////////////
///   MAIN   ///
////////////////


//Initiate Variables
input_path = getDirectory("folder with input images");
roi_path = getDirectory("folder where ROIs are saved"); //cell ROIs
file_list = getFileList(input_path);
ROI_path2 = getDirectory("folder where ROIs for TM vesicles are saved"); //Rab7 ROIs
ROI_path2_files = getFileList(ROI_path2);


for(P = 0; P < file_list.length; P++) {print(file_list[P]);
	if (endsWith(file_list[P], ".czi")) {
	
		run("Bio-Formats Importer", "open=[" + input_path + file_list[P] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension;
		title_raw = getTitle();
		title_wo_MIP = substring(title, 0, lengthOf(title)-4);
		load_ROIs(title_wo_MIP, roi_path);
		cell_no = make_Cell_pics(title_raw);
		for (i = 0; i < cell_no; i++) {
			trackmate(title + "_cell" + i); //this is the image variable
		}
		
		run("Close All");
		roiManager("reset");

	}
}

//////////////////////
///   FUNCTIONS   ///
/////////////////////

function file_in_folder(file_list_tt, file_name) { 
	// checks whether file_name is in the file_list_tt (in case that cell did not have any Rab7-vesicles)
	f_i_f = "False";
	for (k = 0; k < file_list_tt.length; k++) {
		if (file_list_tt[k] == file_name) {
			f_i_f = "True";
			break;
		}
	}
	// returns "True" as string if file is in the list
	// otherwise returns "False"
	return f_i_f;
}


function make_Cell_pics(sel_wind) { 
// function to make duplicates for each cell
	cells = roiManager("Count");
	selectWindow(sel_wind);
	for (j = 0; j < cells; j++){
		selectWindow(sel_wind);
		roiManager("Select", j);
		run("Duplicate...", "title=[" + title + "_cell"+ j + "] duplicate channels=1");
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

	fif = file_in_folder(ROI_path2_files, image + ".zip"); 
	if (fif == "False") {
		// if no ROI file in the folder
		print("For " + image + ": not found");
	}
	
	else {
		// your code when the ROIs are in ROI_path_2
		measureROI(image);
	}        
     
}


function measureROI(image) {
	load_ROIs(image, ROI_path2);
	ROI_amount = roiManager("count");
	threshs = newArray(1500, 1800, 2000, 2300, 2800, 3000); // Array with thresholds to select the best one later

		for (k = 0; k < ROI_amount; k++) {
			
				selectWindow(image + "_sub-bgr");
				roiManager("select", k );
				ROInr = roiManager("index");
				roiManager("rename", "TM_" + ROInr);
				run("Measure");
				
				for (thr = 0; thr < threshs.length; thr++) {
				selectWindow(image + "_sub-bgr");
				roiManager("select", k );
				roiManager("rename", "TM_" + ROInr + "_th_"+ threshs[thr]);
				run("Duplicate...", "title=[" + image + "_ROI-TM" + k + "]");
				run("Clear Outside");
				run("Duplicate...", "title=[" + image + "_ROI-TM_bin" + k + "_th_"+ threshs[thr] + "]");
				run("Threshold...");
				setAutoThreshold("Default dark no-reset");
				setThreshold(threshs[thr], 65535);
				run("Convert to Mask");
				run("Analyze Particles...", "display summarize");			
		}
	}	
			

	roiManager("reset");	
}	
	


function load_ROIs(title, path_roi) {
	roiManager("Open", path_roi + title + ".zip");
}


		
