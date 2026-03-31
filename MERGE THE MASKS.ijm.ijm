// =============================================================
//  PHASE 2 — Batch: merge raw stacks with drawn masks
// =============================================================
macro "Phase 2 - Merge Raw Stacks with Masks" {

    rawDir    = getDirectory("Select 01_RAW/FAST folder");
    maskDir   = getDirectory("Select 03_MASKS folder");
    outputDir = getDirectory("Select 04_MASK_MERGE_AND_BGROUNDSUBTRACT folder");

    list = getFileList(rawDir);
    setBatchMode(true);

    for (i = 0; i < list.length; i++) {
        fname = list[i];

        // Accept .nd2 or .tif, skip any already-merged files
        if (!endsWith(fname, ".nd2") && !endsWith(fname, ".tif") && !endsWith(fname, ".tiff")) continue;
        if (indexOf(fname, "Merged") != -1) continue;

        baseName   = replace(fname, ".nd2", "");
        baseName   = replace(baseName, ".tiff", "");
        baseName   = replace(baseName, ".tif", "");
        maskPath   = maskDir   + baseName + "_masks.tif";
        outputPath = outputDir + "Merged_" + baseName + ".tif";

        if (!File.exists(maskPath)) {
            print("SKIPPED (no mask found): " + baseName);
            continue;
        }
        if (File.exists(outputPath)) {
            print("SKIPPING (output exists): " + baseName);
            continue;
        }

        // --- 1. Open raw stack via Bio-Formats ---
        run("Bio-Formats Importer",
            "open=[" + rawDir + fname + "] color_mode=Default view=Hyperstack stack_order=XYCZT");
        rawID    = getImageID();
        rawTitle = getTitle();
        getDimensions(w, h, c_dim, s_dim, f_dim);
        totalFrames = maxOf(s_dim, f_dim);

        // Read frame interval from embedded metadata (seconds)
        frameInterval = Stack.getFrameInterval();

        //photobleach correction
        run("Bleach Correction", "correction=[Exponential Fit]");

        // Subtract background
        run("Subtract Background...", "rolling=20 stack");

        // Inherit spatial calibration
        getVoxelSize(vw, vh, vd, vunit);

        // --- 2. Open 2-channel mask and split ---
        open(maskPath);
        maskTitle = getTitle();   // e.g. "baseName_masks.tif"
        run("Split Channels");
        // Split produces: "C1-<maskTitle>" (Nuc) and "C2-<maskTitle>" (Cyto)
        nucSingleTitle  = "C1-" + maskTitle;
        cytoSingleTitle = "C2-" + maskTitle;

        // --- 3. Expand each 2-D mask to a full-depth stack ---
        expandToStack(nucSingleTitle,  totalFrames, "NucMask");
        expandToStack(cytoSingleTitle, totalFrames, "CytoMask");

        // --- 4. Build edge stacks BEFORE merge (merge consumes its inputs) ---
        selectWindow("CytoMask");
        run("Duplicate...", "duplicate");
        rename("CytoEdges");
        selectWindow("CytoEdges");
        run("Find Edges", "stack");
        run("Cyan");

        selectWindow("NucMask");
        run("Duplicate...", "duplicate");
        rename("NucEdges");
        selectWindow("NucEdges");
        run("Find Edges", "stack");
        run("Green");

        // --- 5. Merge channels ---
        run("Merge Channels...",
            "c1=[" + rawTitle + "] c2=[NucMask] c3=[CytoMask] create");
        mergedID = getImageID();

        // --- 6. Stamp edge overlays ---
        selectImage(mergedID);
        run("Add Image...", "image=CytoEdges x=0 y=0 opacity=100 zero");
        run("Add Image...", "image=NucEdges  x=0 y=0 opacity=100 zero");

        close("CytoEdges");
        close("NucEdges");

        // --- 7. Set stack properties ---
        selectImage(mergedID);
        Stack.setActiveChannels("100");
        run("Properties...",
            "channels=3 slices=1 frames=" + totalFrames +
            " unit=" + vunit +
            " pixel_width=" + vw +
            " pixel_height=" + vh +
            " voxel_depth=1.0000" +
            " frame=[" + frameInterval + " sec]");

        // --- 8. Save ---
        saveAs("Tiff", outputPath);
        run("Close All");
        print("Done: " + baseName);
    }

    setBatchMode(false);
    print("=== Phase 2 complete ===");
}


// =============================================================
//  Helper: expand a single-frame image into an nFrames stack.
//  Closes the source window; leaves newName as the active image.
// =============================================================
function expandToStack(sourceTitle, nFrames, newName) {
    selectWindow(sourceTitle);
    w = getWidth();
    h = getHeight();
    sourceID = getImageID();

    // Create a blank stack at the target depth
    newImage(newName, "16-bit black", w, h, nFrames);
    stackID = getImageID();

    // Copy the single frame into every slice of the new stack
    for (t = 1; t <= nFrames; t++) {
        selectImage(sourceID);
        run("Select All");
        run("Copy");
        selectImage(stackID);
        setSlice(t);
        run("Paste");
    }
    run("Select None");

    // Close the (now redundant) single-frame source
    selectImage(sourceID);
    close();

    selectWindow(newName);
}