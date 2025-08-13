//E. Gutierrez-Galindo et al, University of Stuttgart (2024)

//This macro was used to segment Rab7/TGN46 vesicles and measure PKD3en-mGFP intensity (MFI) in those ROIs, across all Z-sections of the cell (colocalization)
//Requires: 1) Z-stacks ("input_path"); 2) Cell ROIs (cell outlines) saved in "roi_path" - this can be saved manually or using Save_ROIs macro


function make_Cell_pics(sel_wind) { 
	// function to make duplicates for each cell
	cells = roiManager("Count");
	selectWindow(sel_wind);
	for (j = 0; j < cells; j++){
		selectWindow(sel_wind);
		roiManager("Select", j);
		run("Duplicate...", "title=[" + title + "_cell"+ j + "] duplicate channels=1-3");
		run("Clear Outside", "stack");
	}
	return cells;
}

function analyze_vesicles(image, vesicle_channel, thresh_vesicle, vesicle_name) {
	roiManager("reset");
	selectWindow(image);
	run("Duplicate...", "title=[" + image + "_coloc] duplicate channels=1-3");
	run("Split Channels");

	selectWindow("C3-" + image + "_coloc"); // PKD3-GFP
	run("Duplicate...", "title=[" + image + "_PKD3" + "] duplicate");
	close("C3-" + image + "_coloc");

	selectWindow("C" + vesicle_channel + "-" + image + "_coloc");
	run("Duplicate...", "title=[" + image + "_" + vesicle_name + "] duplicate");
	close("C" + vesicle_channel + "-" + image + "_coloc");
	
	selectWindow(image + "_" + vesicle_name);
	Stack.getDimensions(width, height, channels, slices, frames);

	run("Set Measurements...", "mean integrated area display redirect=None decimal=3");

	selectWindow(image + "_" + vesicle_name);
	Stack.getDimensions(w, h, ch, sl, fr);
	print("Slices detectados en " + image + "___" + vesicle_name + ": " + sl);
	
	// Segment vesicles (TGN46 or Rab7)
	setThreshold(thresh_vesicle, 65535);
	setOption("BlackBackground", false);
	run("Convert to Mask" , "create");
	rename(image + "_" + vesicle_name + "_bin");
	
	
	for (s = 1; s <= slices; s++) {
		print("Slice " + s);

		selectWindow(image + "_" + vesicle_name + "_bin");
		setSlice(s);
		
		run("Duplicate...", "title=" + image + "_bin_" + vesicle_name + "_slice" + s + " use"); //duplicates only stack in use
		selectWindow(image + "_bin_" + vesicle_name + "_slice" + s);
		run("Analyze Particles...", "add slice");

		roi_count = roiManager("count");
		print("ROIs detectadas en " + image + "___" + vesicle_name + "_ slice" + s + ": " + roi_count);

		if (roi_count > 0) {

			// Binarize PKD3-GFP
			selectWindow(image + "_PKD3");
			setSlice(s);
			run("16-bit");
			
			run("Duplicate...", "use"); //duplicates only stack in use
			
			setThreshold(600, 65535); //threshold PKD3-GFP
			run("Convert to Mask");
			rename(image + "bin_PKD3_slice" + s);
			run("Analyze Particles...", "display"); //total segmented PKD3 int density / area

			
			// Colocalization PKD3 ∩ Rab7/Vps35
			selectWindow(image + "bin_PKD3_slice" + s);
			imageCalculator("AND create", image + "bin_PKD3_slice" + s, image + "_bin_" + vesicle_name + "_slice" + s);
			rename(image + "PKD3_coloc" + "_" + vesicle_name + "_slice_" + s);
			run("Analyze Particles...", "display"); //colocalising PKD3 int density/ area
			

			roiManager("reset");
		}
		
		close(image + "bin_PKD3_slice" + s);
		close(image + "_bin_" + vesicle_name + "_slice" + s);
		close(image + "PKD3_coloc" + "_" + vesicle_name + "_slice_" + s);
	}
	
}


function load_ROIs(title, path_roi) {
	roiManager("reset");
	roiManager("Open", path_roi + title + ".zip");
}

/////////////////
///   MAIN   ///
////////////////
input_path = getDirectory("folder with input images");
roi_path = getDirectory("folder where ROIs are saved");
file_list = getFileList(input_path);
output_path = getDirectory("folder where analysed images are saved");

thresh_TGN46 = 7000;
thresh_Rab7 = 8000;

for (P = 0; P < file_list.length; P++) {
	print(file_list[P]);
	if (endsWith(file_list[P], ".czi")) {
		run("Bio-Formats Importer", "open=[" + input_path + file_list[P] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension;
		title_raw = getTitle();
		load_ROIs(title, roi_path);
		cell_no = make_Cell_pics(title_raw);
		for (i = 0; i < cell_no; i++) {
			image_id = title + "_cell" + i;
			analyze_vesicles(image_id, 1, thresh_TGN46, "TGN46");  
			analyze_vesicles(image_id, 2, thresh_Rab7, "Rab7");  
		}
		run("Close All");
		roiManager("reset");
	}
}
