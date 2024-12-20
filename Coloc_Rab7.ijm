//E. Gutierrez-Galindo et al, University of Stuttgart (2024)

//This macro was used to segment Lamp1/CD63/Vps35+ vesicles and measure Rab7 intensity (MFI) in those ROIs, across all Z-sections of the cell (colocalization)
//Requires: 1) Z-stacks ("input_path"); 2) Cell ROIs (cell outlines) saved in "roi_path" - this can be saved manually or using Save_ROIs macro


/////////////////
///   MAIN   ///
////////////////
		
//Initiate Variables
input_path = getDirectory("folder with input images");
roi_path = getDirectory("folder where ROIs are saved");
file_list = getFileList(input_path);
output_path = getDirectory("folder where analysed images are saved")

//Threshold - adjust to signal
thresh_a = 2000; //Lamp1, CD63, Vps35

for(P = 0; P < file_list.length; P++) {print(file_list[P]);
	if (endsWith(file_list[P], ".czi")) {
	
		run("Bio-Formats Importer", "open=[" + input_path + file_list[P] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension;
		title_raw = getTitle();
		load_ROIs(title, roi_path);
		cell_no = make_Cell_pics(title_raw);
		for (i = 0; i < cell_no; i++) {
			coloc_Rab7_int(title + "_cell" + i);
		}
	
		
		run("Close All");
		roiManager("reset");

	}
}


//////////////////////
///   FUNCTIONS   ///
/////////////////////


function make_Cell_pics(sel_wind) { 
// function to make duplicates for each individual cell
	cells = roiManager("Count");
	selectWindow(sel_wind);
	for (j = 0; j < cells; j++){
		selectWindow(sel_wind);
		roiManager("Select", j);
		run("Duplicate...", "title=[" + title + "_cell"+ j + "] duplicate channels=1-2");
		run("Clear Outside", "stack");
	}
	return cells;
}


function coloc_Rab7_int(image) {

//for every Z, analyse Lamp1/CD63/Vps35 vesicles, add to ROI, measure Rab7 MFI in corresponding Z
//Adjust threshold (Thresh_a) above to segment the vesicles
//For each cell, it stops to check if the threshold needs to be corrected to match thresh_a before measuring (sometimes it does not work) - if not correct, adjust manually to thresh_a and then press ok.
	
	selectWindow(image);
	run("Duplicate...", "title=[" + image + "_coloc] duplicate channels=1-2");
//	run("8-bit");
	run("Split Channels");
	selectWindow("C1-" + image + "_coloc");
	run("Duplicate...", "title=[" + image + "_Rab7" + "] duplicate");
	
	
	selectWindow("C2-" + image + "_coloc");
	run("Duplicate...", "title=[" + image + "_Lamp1" + "] duplicate");
	run("Threshold...");
	setThreshold(thresh_a, 65535, "raw");
	setOption("BlackBackground", false);
	
	
	waitForUser("Correct threshold?");
	run("Convert to Mask", "stack create");
	rename(image + "_Lamp1-bin");
	
//	Create output image for Lamp1-bin
	run("Duplicate...", "title=[" + image + "_Lamp1_bin" + "] duplicate");
	saveAs("Tiff", output_path + image + "_Lamp1_bin" + ".tif");

	
	roiManager("reset");
	selectWindow(image + "_Lamp1-bin");
	
    run("Set Measurements...", "area mean integrated display redirect=None decimal=3");
	selectWindow(image + "_Lamp1-bin");
	Stack.getDimensions(width,height,channels,slices,frames);

	for (s = 1; s <= slices; s++) {
		selectWindow(image + "_Lamp1-bin");
		Stack.setSlice(s);		
		ROI_amount = roiManager("count");
    	print(image + "_slice_" + s);
			
		run("Select None");
		run("Analyze Particles...", "add slice");
		ROI_amount_2 = roiManager("count");
		print("ROI_amount_2=" + ROI_amount_2);
			
		if (ROI_amount == ROI_amount_2) {
			print("no coloc found for slice " + image + "_slice_" + s);
			continue;
		}
		
		for (k = ROI_amount; k < ROI_amount_2; k++) {
			roiManager("select", k);
			roiManager("rename", "slice" + s + "_roi" + k);
			selectWindow(image + "_Rab7");
			roiManager("select", k);
			roiManager("measure");
		}
				
	}
	
}
	

function load_ROIs(title, path_roi) {
	roiManager("Open", path_roi + title + ".zip");
}


