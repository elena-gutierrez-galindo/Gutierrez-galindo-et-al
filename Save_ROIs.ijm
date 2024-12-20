//E. Gutierrez-Galindo et al, University of Stuttgart (2024)

//This macro was used to save ROIs of individual cells used for further analysis
//Requires: 1) MIPs or Z-stacks

//1) Draw cell outline using "Freehand selection" (it is recommended to increase Brightness & Contrast)
//2) Add to ROI Manager (T)
//3) Run macro


//path where ROIs should be saved
path = //add path here;

		title = File.nameWithoutExtension;
		title_raw = getTitle();
		title_wo_MIP = substring(title, 0, lengthOf(title)-4);
		print(title_wo_MIP);


roiManager("Save", path + title_wo_MIP + ".zip");
roiManager("reset");

selectWindow(title_raw);
run("Close");
