//E. Gutierrez-Galindo et al, University of Stuttgart (2024)

//This macro was used to measure mean fluorescence intensity (MFI) of each channel, using maximum intensity projections (MIPs)
//Requires: 1) MIPs ("input_path"); 2) Cell ROIs (cell outlines) saved in "roi_path" - this can be saved manually or using Save_ROIs macro


/////////////////
///   MAIN   ///
////////////////


//Initiate Variables
input_path = getDirectory("folder with input images");
roi_path = getDirectory("folder where ROIs are saved");
file_list = getFileList(input_path);


for(P = 0; P < file_list.length; P++) {print(file_list[P]);
	if (endsWith(file_list[P], ".czi")) {
	
		run("Bio-Formats Importer", "open=[" + input_path + file_list[P] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension;
		title_raw = getTitle();
		title_wo_MIP = substring(title, 0, lengthOf(title)-4);
		load_ROIs(title_wo_MIP, roi_path);
		cell_no = make_Cell_pics(title_raw);

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
		//first select the j and then the slice, otherwise it measures the slice used for creating the ROI;
		//use the command run(measure) - not roimanager(measure)
		//omit slices not needed by using //
		
		run("Set Measurements...", "area mean integrated display redirect=None decimal=3");
		selectWindow(sel_wind);
		roiManager("Select", j);
		setSlice(1);
		run("Measure");
		
		roiManager("Select", j);
		setSlice(2);
		run("Measure");
		
		roiManager("Select", j);
		setSlice(3);
		run("Measure");
	}
	roiManager("reset");
	return cells;

}


function load_ROIs(title, path_roi) {
	roiManager("Open", path_roi + title + ".zip");
}

