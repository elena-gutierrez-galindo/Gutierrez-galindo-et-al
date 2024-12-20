//E. Gutierrez-Galindo et al, University of Stuttgart (2024)

//This macro was used to mmeasure colocalization between phospho-PKD (pPKD) and PKD3-GFP vesicles using maximum intensity projections (MIP).
//It creates binary images of both channels and uses the function "Colocalization threshold" to cerate a binary image of the colocalizing vesicles/area; this area is then measured (area, count, intensity)
//Requires: 1) MIPs ("input_path"); 2) Cell ROIs (cell outlines) saved in "roi_path" - this can be saved manually or using Save_ROIs macro


/////////////////
///   MAIN   ///
////////////////

//Initiate Variables
input_path = getDirectory("folder with input images");
roi_path = getDirectory("folder where ROIs are saved");
file_list = getFileList(input_path);
output_path = getDirectory("folder where analysed images are saved")


for(P = 0; P < file_list.length; P++) {print(file_list[P]);
	if (endsWith(file_list[P], ".czi")) {
	
		run("Bio-Formats Importer", "open=[" + input_path + file_list[P] + "] autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		title = File.nameWithoutExtension;
		title_raw = getTitle();
		title_wo_MIP = substring(title, 0, lengthOf(title)-4);
		load_ROIs(title_wo_MIP, roi_path);
		cell_no = make_Cell_pics(title_raw);
		for (i = 0; i < cell_no; i++) {
			pPKD_GFP(title + "_cell" + i);
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
		roiManager("measure");
		run("Duplicate...", "title=[" + title + "_cell"+ j + "] duplicate channels=1-2");
		run("Clear Outside", "stack");
	}
	roiManager("reset");
	return cells;

}

function pPKD_GFP(image) {
		
	//split stack and rename
	negative_cells = 0; //this gives a count of cells where no PKD3 vesicles were found (no active PKD3)
	run("Set Measurements...", "area mean integrated display redirect=None decimal=3");
	selectWindow(image);
	run("8-bit");
	run("Split Channels");
	
	//rename pPKD Channel
	selectWindow("C1-"+image);
	run("Duplicate...", "dup_C1-"+image);
	selectWindow("C1-" + image + "-1");
	rename(image+"_phospho");
	
	//rename GFP channel
	selectWindow("C2-"+image);
	run("Duplicate...", "dup_C2-"+image);
	selectWindow("C2-" + image + "-1");
	rename(image+"_GFP");
	
	//start detecting GFP
		selectWindow(image+"_GFP");
		run("Duplicate...", " ");
		rename(image + "_GFP_coloc");
		selectWindow(image+"_GFP");
		setAutoThreshold("Default dark");
		run("Threshold...");
		setThreshold(15, 255); //this might need to be adapted for each replicate
		setOption("BlackBackground", false);
		run("Convert to Mask");
		rename(image + "_GFP_bin");
		
		//detect phospho mask
		selectWindow(image+"_phospho");
		run("Duplicate...", " ");
		rename(image + "_phospho_coloc");
		selectWindow(image+"_phospho");
		run("Threshold...");
		setThreshold(47, 255); //this might need to be adapted for each replicate
		setOption("BlackBackground", false);
		run("Convert to Mask");
		rename(image + "_phospho_bin");
	
	
		//segment colocalizing area
		run("Colocalization Threshold", "channel_1=["+image+"_phospho_bin] channel_2=["+image+"_GFP_bin] use=[Channel 1] channel=[Red : Green] show include");
		rename(image + "_coloc");
		run("Duplicate...", " ");
				
		saveAs("Tiff", output_path + image+"_coloc" + ".tif");
		
		//measure colocalizing area
		selectWindow(image+"_coloc");
		run("16-bit");
		setThreshold(86, 255);
		rename(image+"_coloc_bin");

		
		ROI_amount = roiManager("count");
		
		run("Analyze Particles...", "pixel show=Outlines display add");
		ROI_amount_2 = roiManager("count");		
		
		if (ROI_amount == ROI_amount_2) {
			print("no coloc found for image " + image);
			negative_cells++;
			continue;
		}
		
				
		
		selectWindow(image + "_GFP_coloc");
		for (k = ROI_amount; k < ROI_amount_2; k++) {
			roiManager("select", k);
			roiManager("measure");	
		}
	

		selectWindow(image + "_phospho_coloc");
		for (k = ROI_amount; k < ROI_amount_2; k++) {
			roiManager("select", k);
			roiManager("measure");	
		}

	roiManager("reset");
	

	}


function load_ROIs(title, path_roi) {
	roiManager("Open", path_roi + title + ".zip");
}

		
