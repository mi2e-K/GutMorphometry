// ============================================================
//  GutMorphometry_.ijm
//  Fiji/ImageJ macro — mouse gut Swiss roll morphometry
//  HE staining, 4x magnification
//
//  Filename convention : AnimalID_Region_Section.tif
//  Region codes        : U (upper SI), M (mid SI),
//                        L (lower SI), C (colon)
//
//  Measurement workflow (structure types):
//    1. Draw line / polyline on the image
//    2. Press T  →  line is registered (stays visible as overlay)
//    3. Repeat for all measurements of this type
//    4. Click OK  →  all newly registered lines are measured at once
//
//  Color coding per measurement type:
//    villus_height  green
//    villus_width   orange
//    crypt_depth    yellow
//    crypt_width    cyan
//    muscularis     magenta
//
//  Output: one CSV per image  (e.g. 622_U_1_results.csv)
//          Data is always APPENDED — existing rows are never
//          overwritten.
// ============================================================

var g_orig_id    = 0;
var g_csv        = "";
var g_json       = "";   // JSON Lines file for line/polygon coordinates
var g_animal     = "";
var g_region     = "";
var g_section    = "";
var g_loop       = "";
var g_script_dir = "";   // directory that contains musc_thickness.py
var g_python_cmd = "python";  // change to full path if needed
                               // e.g. "C:\\Users\\ksmhp\\anaconda3\\python.exe"
var g_scale_file      = "";   // path to saved pixel_scale.txt (set at startup)
var g_do_calibration  = true; // persistent calibration preference (survives loop iterations)
var g_villus_areas    = "";   // deferred villus_area data: "lp|vn|area_um2;" per entry

// ---------------------------------------------------------------
macro "Gut Morphometry [G]" {

    // ----------------------------------------------------------
    // 0. Open image if none is loaded
    // ----------------------------------------------------------
    if (nImages == 0) {
        path = File.openDialog("Select a gut section image (TIF)");
        if (path == "") exit("No image selected. Macro aborted.");
        open(path);
    }

    g_orig_id = getImageID();
    title     = getTitle();

    // ----------------------------------------------------------
    // 1. Parse filename  (e.g. 622_U_1.tif)
    // ----------------------------------------------------------
    dot = lastIndexOf(title, ".");
    if (dot > 0) base = substring(title, 0, dot);
    else         base = title;
    parts = split(base, "_");

    if (parts.length >= 2) {
        g_animal  = parts[0];
        g_region  = toUpperCase(parts[1]);
        if (parts.length >= 3) g_section = parts[2];
        else                   g_section = "1";
    } else {
        Dialog.create("Image Metadata");
        Dialog.addString("Animal ID:", "");
        Dialog.addChoice("Region:", newArray("U","M","L","C"), "U");
        Dialog.addString("Section #:", "1");
        Dialog.show();
        g_animal  = Dialog.getString();
        g_region  = Dialog.getChoice();
        g_section = Dialog.getString();
    }

    // ----------------------------------------------------------
    // 2. Initialise per-image output CSV
    //    File name mirrors the image: 622_U_1_results.csv
    //    Header is written only when the file is new.
    //    All subsequent runs append — existing data is preserved.
    // ----------------------------------------------------------
    img_dir = getDirectory("image");
    out_dir = img_dir + base + File.separator;
    if (!File.exists(out_dir)) File.makeDirectory(out_dir);

    g_csv   = out_dir + base + "_results.csv";

    // Locate musc_thickness.py in Fiji's AutoRun scripts folder.
    // getDirectory("imagej") always returns the Fiji installation root,
    // regardless of where this macro file is stored.
    fiji_dir     = getDirectory("imagej");
    g_script_dir = fiji_dir + "scripts" + File.separator +
                   "Plugins"            + File.separator +
                   "AutoRun"            + File.separator;
    config_dir   = fiji_dir + "GutMorphometry_config" + File.separator;
    if (!File.exists(config_dir)) File.makeDirectory(config_dir);
    g_scale_file = config_dir + "pixel_scale.txt";

    // ---- Clean up files/dirs accidentally placed inside AutoRun --------
    // Fiji tries to execute everything in AutoRun at startup — non-script
    // files (e.g. pixel_scale.txt, musc_params.json) cause spurious Log
    // output.  Migrate any legacy files and delete the originals so the
    // next Fiji start is clean.
    _old_scale = g_script_dir + "pixel_scale.txt";
    if (File.exists(_old_scale)) {
        // Migrate calibration value to new location (don't overwrite if
        // the user has already re-calibrated in the new location).
        if (!File.exists(g_scale_file)) {
            File.saveString(File.openAsString(_old_scale), g_scale_file);
        }
        File.delete(_old_scale);
    }
    _wrong_json = g_script_dir + "musc_params.json";
    _wrong_json2 = g_script_dir + "GutMorphometry_config" + File.separator + "musc_params.json";
    _wrong_dir   = g_script_dir + "GutMorphometry_config";
    if (File.exists(_wrong_json))  File.delete(_wrong_json);
    if (File.exists(_wrong_json2)) File.delete(_wrong_json2);
    if (File.exists(_wrong_dir))   File.delete(_wrong_dir);
    // --------------------------------------------------------------------

    if (!File.exists(g_csv)) {
        File.saveString(
            "AnimalID,Region,Section,LoopPosition,MeasurementType,Value_um,Notes\n",
            g_csv);
    }

    // JSON Lines file for coordinate data (one JSON object per line)
    g_json         = out_dir + base + "_coords.jsonl";
    g_villus_areas = "";

    // ----------------------------------------------------------
    // 3. Session setup  (scale calibration is embedded here)
    // ----------------------------------------------------------
    if (File.exists(g_scale_file)) {
        saved_um_str = File.openAsString(g_scale_file);
        saved_um_px  = parseFloat(saved_um_str);
        scale_opts   = newArray(
            "Use saved  (" + d2s(saved_um_px, 4) + " um/px)",
            "Auto-detect from scale bar ROI",
            "Draw line on scale bar",
            "Enter value (um/px)",
            "Skip — scale already set");
        scale_def = scale_opts[0];
    } else {
        scale_opts = newArray(
            "Auto-detect from scale bar ROI",
            "Draw line on scale bar",
            "Enter value (um/px)",
            "Skip — scale already set");
        scale_def = scale_opts[0];
    }

    Dialog.create("Session Setup — " + g_animal + "_" + g_region + "_" + g_section);
    Dialog.addMessage("Output: " + g_csv);
    Dialog.addChoice("Scale (500 um bar):", scale_opts, scale_def);
    Dialog.addChoice("Loop position in Swiss roll:",
                     newArray("outer","middle","inner"), "outer");
    Dialog.show();

    scale_method = Dialog.getChoice();
    g_loop       = Dialog.getChoice();

    // ----------------------------------------------------------
    // 4. Execute scale calibration based on chosen method
    // ----------------------------------------------------------
    if (startsWith(scale_method, "Use saved")) {
        run("Set Scale...",
            "distance=" + d2s(500.0 / saved_um_px, 6) + " known=500 unit=um global");
        showStatus("Scale: " + d2s(saved_um_px, 4) + " um/px (saved)");

    } else if (scale_method == "Enter value (um/px)") {
        Dialog.create("Enter pixel size");
        Dialog.addNumber("Pixel size (um/px):", 1.0);
        Dialog.show();
        entered_um = Dialog.getNumber();
        run("Set Scale...",
            "distance=" + d2s(500.0 / entered_um, 6) + " known=500 unit=um global");
        File.saveString(d2s(entered_um, 8), g_scale_file);
        showStatus("Scale: " + d2s(entered_um, 4) + " um/px");

    } else if (scale_method == "Auto-detect from scale bar ROI") {
        detected_px = detectScaleBarPx();
        if (detected_px > 0) {
            run("Set Scale...",
                "distance=" + detected_px + " known=500 unit=um global");
            det_um = 500.0 / detected_px;
            File.saveString(d2s(det_um, 8), g_scale_file);
            showStatus("Scale: " + d2s(det_um, 4) + " um/px (auto-detected)");
        }

    } else if (scale_method == "Draw line on scale bar") {
        setTool("line");
        waitForUser("Scale Calibration",
            "Draw a straight line along the full length of\n" +
            "the 500 um scale bar, then click OK.");
        getLine(x1, y1, x2, y2, lw);
        if (x1 == -1) {
            showMessage("Warning", "No line detected — scale was NOT changed.");
        } else {
            line_px = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
            run("Set Scale...",
                "distance=" + line_px + " known=500 unit=um global");
            line_um = 500.0 / line_px;
            File.saveString(d2s(line_um, 8), g_scale_file);
            showStatus("Scale: " + d2s(line_um, 4) + " um/px");
        }
    }
    // "Skip" → do nothing, scale remains as-is

    // ----------------------------------------------------------
    // 5. Open ROI Manager
    //    Lines registered with T remain visible as overlays
    //    throughout the session, building a visual record.
    // ----------------------------------------------------------
    if (!isOpen("ROI Manager")) run("ROI Manager...");
    roiManager("show all");

    // ----------------------------------------------------------
    // 6. Measurement menu by region
    // ----------------------------------------------------------
    is_small = (startsWith(g_region, "U") || startsWith(g_region, "M") || startsWith(g_region, "L"));

    if (is_small) {
        mtype_list = newArray(
            "villus",
            "crypt_depth",
            "muscularis",
            "--- Save and Exit ---"
        );
    } else {
        mtype_list = newArray(
            "crypt_depth",
            "muscularis",
            "--- Save and Exit ---"
        );
    }

    // ----------------------------------------------------------
    // 7. Main measurement loop
    // ----------------------------------------------------------
    while (true) {
        Dialog.create("Gut Morphometry — " + g_animal + " | " + g_region);
        Dialog.addChoice("Measurement type:", mtype_list, mtype_list[0]);
        Dialog.addChoice("Loop position:", newArray("outer","middle","inner"), g_loop);
        Dialog.addNumber("Measurements per session (villus / crypt):", 10);
        Dialog.addCheckbox("Pick background colour (villus width segmentation)", false);
        Dialog.addString("Notes (optional):", "");
        Dialog.addCheckbox("Calibrate muscularis threshold (ignored for other types)",
                           g_do_calibration);
        Dialog.show();

        mtype            = Dialog.getChoice();
        g_loop           = Dialog.getChoice();
        n_session        = maxOf(1, round(Dialog.getNumber()));
        do_bg_pick       = Dialog.getCheckbox();
        notes            = Dialog.getString();
        g_do_calibration = Dialog.getCheckbox();

        if (mtype == "--- Save and Exit ---") {
            saveOverlayImage();
            break;
        }

        if (startsWith(mtype, "muscularis")) {
            measureMuscularis(mtype, notes);
        } else if (mtype == "villus") {
            measureVillus(notes, n_session, do_bg_pick);
        } else {
            measureBatch(mtype, notes, n_session);
        }

        selectImage(g_orig_id);
    }

    // Completion message is shown inside saveOverlayImage().
}

// ============================================================
//  measureBatch
//
//  Batch measurement using the ROI Manager (T-key workflow).
//  After the user presses OK, all newly registered lines are:
//    - Coloured by measurement type
//    - Numbered sequentially (e.g. CD-1, CD-2 ...)
//    - Measured and appended to the CSV
//
//  Previously registered overlays are preserved but not
//  re-measured (tracked via n_before / n_after).
// ============================================================
function measureBatch(mtype, notes, n_target) {

    is_tortuous = endsWith(mtype, "height");   // depth uses straight line

    if (is_tortuous) {
        setTool("polyline");
        draw_msg =
            "Trace along the structure:\n" +
            "  - Click to add waypoints\n" +
            "  - Double-click or right-click to finish";
    } else {
        setTool("line");
        draw_msg = "Draw a straight line across the structure.";
    }

    n_before   = roiManager("count");
    start_num  = countType(mtype) + 1;   // continue numbering from CSV
    color      = getColor(mtype);
    prefix     = getPrefix(mtype);
    n_roi_prev = n_before;
    n_measured = 0;

    // Per-line loop: each iteration = one line.
    // Break when n_target reached or user clicks OK without drawing.
    while (true) {
        if (n_measured >= n_target) break;

        cur_num = start_num + n_measured;

        waitForUser(
            prefix + "-" + cur_num + "  (" + (n_measured + 1) + "/" + n_target + ")",
            draw_msg + "\n\n" +
            "Press T to register the line.\n" +
            "Then click OK.\n\n" +
            "► Click OK WITHOUT drawing to stop early.");

        n_roi_now = roiManager("count");

        if (n_roi_now == n_roi_prev) {
            // No new ROI — check for forgotten T
            sel = selectionType();
            while (sel == 5 || sel == 6 || sel == 7) {
                waitForUser(
                    "Forgot to press T?  (" + (n_measured + 1) + "/" + n_target + ")",
                    "A line is drawn but not yet registered.\n\n" +
                    "► Press T to register it, then click OK.\n\n" +
                    "► Click Cancel to end measurements.");
                if (roiManager("count") > n_roi_prev) break;
                sel = selectionType();
            }
            if (roiManager("count") == n_roi_prev) break;  // truly done
            n_roi_now = roiManager("count");
        }

        // Label, color, measure, and save coords for each newly registered ROI
        for (j = n_roi_prev; j < n_roi_now; j++) {
            line_lbl = prefix + "-" + (start_num + j - n_before);
            roiManager("select", j);
            roiManager("rename", line_lbl);
            Roi.setStrokeColor(color);
            Roi.setStrokeWidth(2);
            roiManager("update");

            getSelectionCoordinates(xc, yc);

            run("Measure");
            val = getResult("Length", nResults - 1);
            IJ.deleteRows(nResults - 1, nResults - 1);

            appendCSV(mtype, val, notes);
            appendJSON(mtype, line_lbl, val, notes, xc, yc);
            n_measured++;
        }

        n_roi_prev = n_roi_now;
        roiManager("show all");
    }

    n_new = n_roi_prev - n_before;
    if (n_new == 0) {
        showMessage("Nothing recorded",
            "No new lines were added to ROI Manager.\n" +
            "Remember to press T after drawing each line.");
    } else {
        roiManager("deselect");
        roiManager("show all");
        showStatus(n_new + " x " + mtype + " recorded.");
    }
    _transferToOverlayAndReset();
}

// ============================================================
//  measureMuscularis
//
//  Automatic thickness via Python (musc_thickness.py / scipy EDT).
//
//  1. User draws a POLYGON tightly around the muscle layer.
//  2. Pixels outside the polygon are set to white (background).
//     The polygon boundary IS the spatial constraint.
//  3. Masked crop → Python → mean EDT thickness (μm).
//  4. Python saves a verification PNG; user confirms before CSV.
//
//  Requires: Python 3 with numpy, scikit-image, scipy, matplotlib.
//            musc_thickness.py in Fiji/scripts/Plugins/AutoRun/
// ============================================================
function measureMuscularis(mtype, notes) {

    // Work on a duplicate so the original image is never modified.
    selectImage(g_orig_id);
    run("Duplicate...", "title=_work_copy");
    work_id = getImageID();

    // ----------------------------------------------------------
    // 1. 500 µm grid overlay — scale reference while drawing
    //    Drawn on the duplicate; removed before cropping so it
    //    never reaches the TIF that Python reads.
    // ----------------------------------------------------------
    selectImage(work_id);
    getVoxelSize(vw, vh, vd, vunit);
    grid_px = round(500.0 / vw);
    iw = getWidth();
    ih = getHeight();
    run("Remove Overlay");
    if (grid_px > 5 && grid_px < iw) {
        for (gx = grid_px; gx < iw; gx += grid_px) {
            makeLine(gx, 0, gx, ih);
            Roi.setStrokeColor("#a600ffff");
            Roi.setStrokeWidth(1);
            Overlay.addSelection();
        }
    }
    if (grid_px > 5 && grid_px < ih) {
        for (gy = grid_px; gy < ih; gy += grid_px) {
            makeLine(0, gy, iw, gy);
            Roi.setStrokeColor("#a600ffff");
            Roi.setStrokeWidth(1);
            Overlay.addSelection();
        }
    }
    run("Select None");

    setTool("polygon");
    waitForUser("Step 1 — Draw ROI around muscle layer",
        "Cyan grid = 500 µm spacing (scale reference).\n" +
        "Cursor µm position is shown live at the bottom of the Fiji window.\n\n" +
        "Draw a POLYGON tightly enclosing the muscle layer.\n" +
        "  - Click to place each vertex\n" +
        "  - Double-click (or right-click) to close the polygon\n\n" +
        "Include a small margin of tissue above and below.\n\n" +
        "Click OK when the polygon is placed.");

    if (selectionType() != 2) {
        showMessage("Skipped",
            "A polygon selection is required — measurement not recorded.");
        if (isOpen(work_id)) { selectImage(work_id); close(); }
        selectImage(g_orig_id);
        return;
    }

    // Remove grid overlay before duplicating — must not be embedded in TIF.
    selectImage(work_id);
    run("Remove Overlay");

    run("Duplicate...", "title=_musc_crop");
    crop_id = getImageID();
    // vw is already set above; re-read on crop to be safe
    getVoxelSize(vw, vh, vd, vunit);

    // Save polygon coordinates for JSON before selection is cleared
    getSelectionCoordinates(mu_xc, mu_yc);

    // Apply polygon mask: set everything outside the ROI to white.
    // Python's tissue-mask (OD threshold) will naturally exclude those
    // white pixels — the polygon boundary becomes the spatial constraint.
    setBackgroundColor(255, 255, 255);
    run("Clear Outside");
    run("Select None");

    if (isOpen(work_id)) { selectImage(work_id); close(); }

    // ----------------------------------------------------------
    // 2. Save masked crop → call Python (thickness + visualization)
    //    Uses cmd.exe to ensure Anaconda Python is called,
    //    not Fiji's bundled Jython.
    // ----------------------------------------------------------
    tmp_img   = getDirectory("temp") + "musc_crop_tmp.tif";
    tmp_out   = getDirectory("temp") + "musc_result_tmp.txt";
    tmp_viz   = getDirectory("temp") + "musc_viz_tmp.png";
    py_script = g_script_dir + "musc_thickness.py";

    if (!File.exists(py_script)) {
        showMessage("Error",
            "Python script not found:\n" + py_script + "\n\n" +
            "Place musc_thickness.py in:\n" + g_script_dir);
        if (isOpen(crop_id)) { selectImage(crop_id); close(); }
        selectImage(g_orig_id);
        return;
    }

    selectImage(crop_id);
    getDimensions(crop_w, crop_h, crop_ch, crop_sl, crop_fr);
    saveAs("Tiff", tmp_img);
    saved_id = getImageID();

    // ----------------------------------------------------------
    // 2b. Threshold calibration (if selected in Step 1)
    // ----------------------------------------------------------
    tune_script = g_script_dir + "tune_muscularis.py";
    if (g_do_calibration) {
        if (File.exists(tune_script)) {
            tune_cmd = g_python_cmd +
                       " \"" + tune_script + "\"" +
                       " \"" + tmp_img     + "\"" +
                       " "   + d2s(vw, 8);
            exec("cmd.exe", "/c", tune_cmd);
        } else {
            showMessage("Warning",
                "tune_muscularis.py not found:\n" + tune_script +
                "\nProceeding without calibration.");
        }
    }

    // ----------------------------------------------------------
    // 2c. Build label and result image path
    //     Label: AnimalID_Region_Section_Loop_MUn
    //     Saved permanently next to the image CSV.
    // ----------------------------------------------------------
    mu_num     = countType(mtype) + 1;
    lbl        = g_animal + "_" + g_region + "_" + g_section +
                 "_" + g_loop + "_MU" + mu_num;
    result_png = out_dir + lbl + "_muscularis.png";

    // cmd.exe /c forces Windows shell resolution → uses PATH Python
    // (avoids Fiji's bundled Jython being called instead)
    py_cmd = g_python_cmd +
             " \"" + py_script  + "\"" +
             " \"" + tmp_img    + "\"" +
             " "   + d2s(vw, 8) +
             " 0.0" +
             " \"" + tmp_out    + "\"" +
             " \"" + tmp_viz    + "\"" +
             " \"" + lbl        + "\"" +
             " \"" + result_png + "\"";
    exec("cmd.exe", "/c", py_cmd);

    File.delete(tmp_img);
    if (isOpen(saved_id)) { selectImage(saved_id); close(); }

    // ----------------------------------------------------------
    // 3. Parse thickness result
    // ----------------------------------------------------------
    if (!File.exists(tmp_out)) {
        showMessage("Error",
            "Python script produced no output.\n" +
            "Ensure these packages are installed:\n" +
            "  pip install numpy scikit-image scipy matplotlib");
        selectImage(g_orig_id);
        return;
    }

    result_str = File.openAsString(tmp_out);
    // Output format: "mean_um,length_um"  (length = area / mean_thickness)
    result_parts = split(result_str, ",");
    mean_um = parseFloat(result_parts[0]);
    if (result_parts.length > 1) {
        length_um = parseFloat(result_parts[1]);
    } else {
        length_um = 0;
    }
    File.delete(tmp_out);

    if (isNaN(mean_um) || mean_um <= 0) {
        showMessage("Warning",
            "Unexpected result from Python: \"" + result_str + "\"\n" +
            "Measurement not recorded.");
        selectImage(g_orig_id);
        return;
    }

    // ----------------------------------------------------------
    // 4. Visual verification — open the 4-panel PNG
    //    Show only when calibration was NOT run
    //    (calibration already provides a real-time preview).
    // ----------------------------------------------------------
    if (!g_do_calibration && File.exists(tmp_viz)) {
        open(tmp_viz);
        viz_id = getImageID();
        waitForUser("Step 3 — Verify detection",
            "Panel 1 : original crop\n" +
            "Panel 2 : detected muscle mask (green overlay)\n" +
            "Panel 3 : local thickness heatmap\n" +
            "Panel 4 : E-H feature map + boundary (yellow)\n\n" +
            "Mean thickness = " + d2s(mean_um, 1) + " um\n" +
            "Measured length = " + d2s(length_um, 0) + " um\n\n" +
            "Result image saved to:\n" + result_png + "\n\n" +
            "Click OK to record this measurement.");
        if (isOpen(viz_id)) { selectImage(viz_id); close(); }
    }
    if (File.exists(tmp_viz)) File.delete(tmp_viz);

    // ----------------------------------------------------------
    // 5. Append to CSV
    // ----------------------------------------------------------
    len_note  = "medial_axis_EDT;eff_len=" + d2s(length_um, 0) + "um";
    if (lengthOf(notes) > 0) {
        auto_note = notes + ";" + len_note;
    } else {
        auto_note = len_note;
    }
    appendCSV(mtype, mean_um, auto_note);
    appendJSON(mtype, lbl, mean_um, auto_note, mu_xc, mu_yc);
    showStatus("Recorded: " + mtype + " = " + d2s(mean_um, 1) + " um  len=" + d2s(length_um, 0) + " um");

    selectImage(g_orig_id);
}

// ============================================================
//  countType
//  Counts how many rows of a given mtype already exist in the
//  CSV, used to continue sequential numbering across sessions.
// ============================================================
function countType(mtype) {
    if (!File.exists(g_csv)) return 0;
    content = File.openAsString(g_csv);
    lines   = split(content, "\n");
    count   = 0;
    for (i = 1; i < lines.length; i++) {   // skip header row
        if (lengthOf(lines[i]) == 0) continue;  // skip empty lines (trailing newline)
        fields = split(lines[i], ",");
        if (fields.length > 4 && fields[4] == mtype) count++;
    }
    return count;
}

// ============================================================
//  getColor / getPrefix
//  Per-type overlay colour and short label prefix.
// ============================================================
function getColor(mtype) {
    if (mtype == "villus_height") return "green";
    if (mtype == "villus_width")  return "orange";
    if (mtype == "crypt_depth")   return "yellow";
    if (mtype == "crypt_width")   return "cyan";
    if (mtype == "muscularis")    return "magenta";
    return "white";
}

function getPrefix(mtype) {
    if (mtype == "villus_height") return "VH";
    if (mtype == "villus_width")  return "VW";
    if (mtype == "crypt_depth")   return "CD";
    if (mtype == "crypt_width")   return "CW";
    if (mtype == "muscularis")    return "MU";
    return "XX";
}

// ============================================================
//  appendCSV
//  Appends one row. Uses global session variables.
//  Commas in Notes are replaced with semicolons.
// ============================================================
function appendCSV(mtype, value, notes) {
    notes = replace(notes, ",", ";");
    row = g_animal + "," + g_region  + "," + g_section + "," +
          g_loop   + "," + mtype     + "," + d2s(value, 2) + "," +
          notes;
    // File.append throws a Java IOException on file-lock errors rather than
    // returning false, so the macro aborts before any error check is reached.
    // Running the write inside a Rhino JavaScript snippet lets us catch the
    // exception and show a helpful message instead of a hard crash.
    csv_js = replace(g_csv, "\\", "\\\\");
    row_js = replace(row,   "\\", "\\\\");
    row_js = replace(row_js, "'",  "\\'");
    js = "try {" +
         "  var fw = new java.io.FileWriter('" + csv_js + "', true);" +
         "  fw.write('" + row_js + "\\n');" +
         "  fw.close();" +
         "  'ok'" +
         "} catch (e) { String(e) }";
    write_result = eval("script", js);
    if (write_result != "ok") {
        showMessage("CSV Write Error",
            "Could not write to:\n  " + g_csv + "\n\n" +
            "The file may be open in Excel or another program.\n" +
            "Please close it and retry.\n\n" +
            "Data that could not be saved:\n  " + row + "\n\n" +
            "Error: " + write_result);
    }
}

// ============================================================
//  measureVillus
//
//  Combined villus height + width measurement (one function, one villus).
//
//  User input (two steps per villus):
//    1. Segmented line along the centerline  →  T to register
//       Arc length = villus_height (measured directly in Fiji)
//       Same polyline = centerline for width calculation in Python
//    2. Polygon ROI tracing the villus outline  →  click OK
//       If cancelled: height is still recorded; width is skipped.
//
//  CSV rows appended:
//    villus_height  value = arc length (µm)
//    villus_width   value = mean_mid width (µm), Notes = w25/50/75/max/n
//
//  QC image: saved as <label>_villus.png next to the CSV.
// ============================================================
function measureVillus(notes, n_target, do_bg) {

    py_script = g_script_dir + "villus_width.py";

    if (do_bg) _pickBgRgb();

    // ----------------------------------------------------------
    // Continuous measurement loop — one iteration = one villus.
    // Break conditions:
    //   1. n_measured reaches n_target (automatic)
    //   2. User clicks OK at Step 1 without drawing (manual early exit)
    // ----------------------------------------------------------
    n_measured = 0;
    while (true) {
        if (n_measured >= n_target) break;

        // Re-read counter each iteration so numbering follows CSV
        vn  = countType("villus_height") + 1;
        lbl = g_animal + "_" + g_region + "_" + g_section +
              "_" + g_loop + "_V" + vn;

        // ----------------------------------------------------------
        // Step 1: centerline (segmented polyline)
        //   arc length  = villus_height
        //   coordinates = centerline for Python width calculation
        // ----------------------------------------------------------
        n_before_cl = roiManager("count");
        setTool("polyline");
        waitForUser(
            "Villus-" + vn + "  Step 1 — Draw centerline  (" + (n_measured + 1) + "/" + n_target + ")",
            "Draw a SEGMENTED LINE along the villus centerline.\n" +
            "  Start at the BASE, end at the TIP.\n\n" +
            "Press T to register.\n" +
            "Then click OK.\n\n" +
            "► Click OK WITHOUT drawing to stop early.");

        if (roiManager("count") == n_before_cl) {
            // No new ROI — check if a line was drawn but T was not pressed.
            // waitForUser keeps image focus so T still works.
            // To end session: press Esc (clears selection) then OK.
            sel = selectionType();
            while (sel == 5 || sel == 6 || sel == 7) {
                // 5=straight line  6=segmented line  7=freehand line
                waitForUser(
                    "Forgot to press T?  (" + (n_measured + 1) + "/" + n_target + ")",
                    "A line is drawn but not yet registered.\n\n" +
                    "► Press T to register it, then click OK.\n\n" +
                    "► To end session: click Cancel.");
                if (roiManager("count") > n_before_cl) break;  // T was pressed
                sel = selectionType();  // Esc clears → sel becomes -1 → loop exits
            }
            if (roiManager("count") == n_before_cl) break;  // truly done
        }

        roi_cl_idx = roiManager("count") - 1;
        roiManager("select", roi_cl_idx);
        roiManager("rename", "VH-" + vn);
        Roi.setStrokeColor("green");
        Roi.setStrokeWidth(2);
        roiManager("update");
        getSelectionCoordinates(cl_x, cl_y);

        // Measure arc length (= villus height)
        run("Measure");
        height_um = getResult("Length", nResults - 1);
        IJ.deleteRows(nResults - 1, nResults - 1);

        // ----------------------------------------------------------
        // Step 2: polygon outline (for width)
        // ----------------------------------------------------------
        setTool("polygon");
        waitForUser(
            "Villus-" + vn + "  Step 2 — Draw outline  (height = " + d2s(height_um, 1) + " µm)",
            "Draw a POLYGON tracing the villus boundary.\n" +
            "Trace close to the tissue edge for accurate width.\n\n" +
            "Click OK when done (no need to press T).\n\n" +
            "Click OK without drawing to record height only and skip width.");

        // Record height regardless of whether polygon was drawn
        h_note = "label=VH-" + vn;
        if (lengthOf(notes) > 0) h_note = notes + ";" + h_note;
        appendCSV("villus_height", height_um, h_note);
        appendJSON("villus_height", "VH-" + vn, height_um, h_note, cl_x, cl_y);
        n_measured++;

        if (selectionType() != 2) {
            // No polygon → skip area and width

            showStatus("VH-" + vn + " height recorded: " + d2s(height_um, 1) + " µm  (width skipped)");
            roiManager("deselect");
            roiManager("show all");
            selectImage(g_orig_id);
            continue;   // move to next villus
        }

        getSelectionCoordinates(poly_x, poly_y);

        // Shoelace formula for polygon area in µm²
        selectImage(g_orig_id);
        getVoxelSize(_vw, _vh, _vd, _vu);
        _area_px2 = 0;
        _np = lengthOf(poly_x);
        for (_k = 0; _k < _np - 1; _k++) {
            _area_px2 += poly_x[_k] * poly_y[_k+1] - poly_x[_k+1] * poly_y[_k];
        }
        _area_px2 += poly_x[_np-1] * poly_y[0] - poly_x[0] * poly_y[_np-1];
        area_um2 = abs(_area_px2) / 2.0 * _vw * _vw;
        // Defer to Save & Exit so all areas appear together, sorted by VH number
        g_villus_areas = g_villus_areas + g_loop + "|" + vn + "|" + d2s(area_um2, 4) + ";";

        // Add polygon to overlay
        roiManager("add");
        roi_poly_idx = roiManager("count") - 1;
        roiManager("select", roi_poly_idx);
        roiManager("rename", "VW-" + vn + "-outline");
        Roi.setStrokeColor("#a600ffff");
        Roi.setStrokeWidth(1);
        roiManager("update");

        if (!File.exists(py_script)) {
            showMessage("Error",
                "villus_width.py not found:\n" + py_script + "\n\n" +
                "Height was recorded; width skipped.");
            selectImage(g_orig_id);
            continue;
        }

        // ----------------------------------------------------------
        // Step 3: tight crop from polygon bounding box
        // ----------------------------------------------------------
        min_x = poly_x[0]; max_x = poly_x[0];
        min_y = poly_y[0]; max_y = poly_y[0];
        for (k = 1; k < lengthOf(poly_x); k++) {
            if (poly_x[k] < min_x) min_x = poly_x[k];
            if (poly_x[k] > max_x) max_x = poly_x[k];
            if (poly_y[k] < min_y) min_y = poly_y[k];
            if (poly_y[k] > max_y) max_y = poly_y[k];
        }
        margin = 30;
        cx = maxOf(0, floor(min_x) - margin);
        cy = maxOf(0, floor(min_y) - margin);
        cw = minOf(getWidth(),  floor(max_x) + 1 + margin) - cx;
        ch = minOf(getHeight(), floor(max_y) + 1 + margin) - cy;

        selectImage(g_orig_id);
        makeRectangle(cx, cy, cw, ch);
        run("Duplicate...", "title=_vw_crop");
        crop_id = getImageID();
        run("Select None");
        getVoxelSize(vw_px, vh_px, vd_px, vunit_px);

        // ----------------------------------------------------------
        // Step 4: save temp files
        // ----------------------------------------------------------
        tmp_img  = getDirectory("temp") + "vw_crop_tmp.tif";
        tmp_poly = getDirectory("temp") + "vw_poly_tmp.json";
        tmp_cl   = getDirectory("temp") + "vw_cl_tmp.json";
        tmp_out  = getDirectory("temp") + "vw_result_tmp.txt";
        qc_png   = out_dir + lbl + "_villus.png";

        selectImage(crop_id);
        saveAs("Tiff", tmp_img);
        saved_id = getImageID();
        selectImage(saved_id); close();

        // Coordinate JSON strings shifted to crop space
        poly_json = "[";
        for (k = 0; k < lengthOf(poly_x); k++) {
            if (k > 0) poly_json += ",";
            poly_json += "[" + d2s(poly_x[k] - cx, 1) + "," + d2s(poly_y[k] - cy, 1) + "]";
        }
        poly_json += "]";

        cl_json = "[";
        for (k = 0; k < lengthOf(cl_x); k++) {
            if (k > 0) cl_json += ",";
            cl_json += "[" + d2s(cl_x[k] - cx, 1) + "," + d2s(cl_y[k] - cy, 1) + "]";
        }
        cl_json += "]";

        poly_js = replace(tmp_poly, "\\", "\\\\");
        cl_js   = replace(tmp_cl,   "\\", "\\\\");
        eval("script",
            "try{var fw=new java.io.FileWriter('" + poly_js + "',false);" +
            "fw.write('" + replace(poly_json,"'","\\'") + "');fw.close();'ok'}" +
            "catch(e){String(e)}");
        eval("script",
            "try{var fw=new java.io.FileWriter('" + cl_js + "',false);" +
            "fw.write('" + replace(cl_json,"'","\\'") + "');fw.close();'ok'}" +
            "catch(e){String(e)}");

        // ----------------------------------------------------------
        // Step 5: run Python (width calculation)
        // ----------------------------------------------------------
        py_cmd = g_python_cmd +
                 " \"" + py_script  + "\"" +
                 " \"" + tmp_img    + "\"" +
                 " "   + d2s(vw_px, 8) +
                 " \"" + tmp_poly   + "\"" +
                 " \"" + tmp_cl     + "\"" +
                 " \"" + tmp_out    + "\"" +
                 " \"" + qc_png     + "\"" +
                 " \"" + lbl        + "\"";
        exec("cmd.exe", "/c", py_cmd);

        File.delete(tmp_img);
        File.delete(tmp_poly);
        File.delete(tmp_cl);

        // ----------------------------------------------------------
        // Step 6: parse results
        // ----------------------------------------------------------
        if (!File.exists(tmp_out)) {
            showMessage("Error",
                "Python script produced no output.\n" +
                "Height was recorded; width skipped.\n\n" +
                "Ensure numpy, scipy, scikit-image, matplotlib are installed.");
            selectImage(g_orig_id);
            continue;
        }

        result_str = File.openAsString(tmp_out);
        rp         = split(result_str, ",");
        File.delete(tmp_out);

        mean_mid = parseFloat(rp[0]);
        w20      = parseFloat(rp[1]);
        w50      = parseFloat(rp[2]);
        w80      = parseFloat(rp[3]);
        wmax     = parseFloat(rp[4]);
        n_valid  = parseFloat(rp[5]);
        qc_flag  = parseFloat(rp[6]);

        if (isNaN(mean_mid) || mean_mid <= 0) {
            showMessage("Warning",
                "Width calculation failed: \"" + result_str + "\"\n" +
                "Height was recorded. Check QC image:\n" + qc_png);
            selectImage(g_orig_id);
            continue;
        }

        // ----------------------------------------------------------
        // Step 7: QC overlay confirmation
        // ----------------------------------------------------------
        if (File.exists(qc_png)) {
            open(qc_png);
            qc_id = getImageID();
            waitForUser("Villus-" + vn + " — Verify width",
                "QC overlay: ROI boundary (cyan) · centerline (yellow)\n" +
                "width chords 20/50/80% (green) · invalid points (red)\n\n" +
                "height   = " + d2s(height_um, 1) + " µm\n" +
                "mean_mid = " + d2s(mean_mid, 1) + " µm\n" +
                "w20=" + d2s(w20, 1) + "  w50=" + d2s(w50, 1) +
                "  w80=" + d2s(w80, 1) + "  wmax=" + d2s(wmax, 1) + "\n" +
                "n_valid=" + d2s(n_valid, 0) + "  qc_flag=" + d2s(qc_flag, 0) + "\n\n" +
                "Click OK to record width and continue to the next villus.");
            if (isOpen(qc_id)) { selectImage(qc_id); close(); }
        }

        // ----------------------------------------------------------
        // Step 8: write width to CSV and JSON
        // ----------------------------------------------------------
        w_note = "w20=" + d2s(w20, 1) + ";w50=" + d2s(w50, 1) +
                 ";w80=" + d2s(w80, 1) + ";wmax=" + d2s(wmax, 1) +
                 ";n=" + d2s(n_valid, 0) + ";pair=VH-" + vn;
        if (lengthOf(notes) > 0) w_note = notes + ";" + w_note;

        appendCSV("villus_width", mean_mid, w_note);
        appendJSON("villus_width", "VW-" + vn, mean_mid, w_note, poly_x, poly_y);

        showStatus("Villus-" + vn + "  height=" + d2s(height_um, 1) +
                   " µm  width_mid=" + d2s(mean_mid, 1) + " µm");
        roiManager("deselect");
        roiManager("show all");
        selectImage(g_orig_id);

    } // end while — villus loop

    roiManager("deselect");
    roiManager("show all");
    selectImage(g_orig_id);
    _transferToOverlayAndReset();
}

// ============================================================
//  _pickBgRgb
//  Lets the user click once on a background region to sample the
//  RGB colour of the slide glass / unstained area.
//  Saves the result as bg_rgb in GutMorphometry_config/musc_params.json
//  so that villus_width.py can use OD-distance tissue segmentation.
// ============================================================
function _pickBgRgb() {
    selectImage(g_orig_id);
    setTool("point");
    waitForUser("Pick background colour",
        "Click on a BACKGROUND area — white/glass region with no tissue.\n" +
        "Then click OK.");

    if (selectionType() != 10) {
        showMessage("Background skipped",
            "No point selected — background colour was not updated.\n" +
            "The existing musc_params.json value (if any) will be used.");
        run("Select None");
        return;
    }

    getSelectionCoordinates(xc, yc);
    pv   = getPixel(round(xc[0]), round(yc[0]));
    bg_r = (pv >> 16) & 0xff;
    bg_g = (pv >>  8) & 0xff;
    bg_b =  pv        & 0xff;
    run("Select None");

    // Save bg_rgb to musc_params.json (update field; preserve other keys)
    config_dir  = getDirectory("imagej") + "GutMorphometry_config" + File.separator;
    params_path = config_dir + "musc_params.json";
    if (!File.exists(config_dir)) File.makeDirectory(config_dir);

    existing = "{}";
    if (File.exists(params_path)) existing = File.openAsString(params_path);

    // Escape for safe embedding inside a JavaScript string literal
    esc = replace(existing, "\\", "\\\\");
    esc = replace(esc, "'",  "\\'");
    esc = replace(esc, "\n", "\\n");
    esc = replace(esc, "\r", "");

    params_js = replace(params_path, "\\", "\\\\");
    js = "try {" +
         "  var p;" +
         "  try { p = JSON.parse('" + esc + "'); } catch(e2) { p = {}; }" +
         "  p.bg_rgb = [" + bg_r + "," + bg_g + "," + bg_b + "];" +
         "  var fw = new java.io.FileWriter('" + params_js + "', false);" +
         "  fw.write(JSON.stringify(p, null, 2));" +
         "  fw.close();" +
         "  'ok'" +
         "} catch(e) { String(e) }";
    result = eval("script", js);

    if (result == "ok") {
        showStatus("Background: R=" + bg_r + " G=" + bg_g + " B=" + bg_b +
                   "  saved to musc_params.json");
    } else {
        showMessage("Warning",
            "Could not save background colour:\n  " + params_path + "\n\n" +
            "Error: " + result);
    }
}

// ============================================================
//  appendJSON
//  Appends one JSON object (one line) to the .jsonl coordinate file.
//  xc / yc are pixel coordinate arrays from getSelectionCoordinates().
//  Both pixel and µm coordinates are stored.
// ============================================================
function appendJSON(mtype, lbl, val_um, notes, xc, yc) {
    if (lengthOf(g_json) == 0) return;

    getVoxelSize(vw, vh, vd, vunit);
    n_pts = lengthOf(xc);

    // Build coordinate arrays as JSON strings
    coords_px = "[";
    coords_um = "[";
    for (k = 0; k < n_pts; k++) {
        if (k > 0) { coords_px += ","; coords_um += ","; }
        coords_px += "[" + d2s(xc[k],      1) + "," + d2s(yc[k],      1) + "]";
        coords_um += "[" + d2s(xc[k] * vw, 2) + "," + d2s(yc[k] * vw, 2) + "]";
    }
    coords_px += "]";
    coords_um += "]";

    notes_safe = replace(notes, "\"", "'");

    json_line = "{" +
        "\"animal\":\"" + g_animal  + "\"," +
        "\"region\":\"" + g_region  + "\"," +
        "\"section\":\"" + g_section + "\"," +
        "\"loop\":\"" + g_loop    + "\"," +
        "\"type\":\"" + mtype     + "\"," +
        "\"label\":\"" + lbl       + "\"," +
        "\"length_um\":" + d2s(val_um, 4) + "," +
        "\"notes\":\"" + notes_safe + "\"," +
        "\"coords_px\":" + coords_px + "," +
        "\"coords_um\":" + coords_um + "}";

    // Use Rhino JavaScript for exception-safe file writing
    json_js = replace(g_json,     "\\", "\\\\");
    line_js = replace(json_line,  "\\", "\\\\");
    line_js = replace(line_js,    "'",  "\\'");
    js = "try {" +
         "  var fw = new java.io.FileWriter('" + json_js + "', true);" +
         "  fw.write('" + line_js + "\\n');" +
         "  fw.close(); 'ok'" +
         "} catch(e) { String(e) }";
    eval("script", js);
}

// ============================================================
//  saveOverlayImage
//  Saves a JPEG of the original image with all ROI overlays
//  burned in.  Called automatically on "Save and Exit".
//  JPEG quality 65% gives a compact confirmation image.
// ============================================================
function saveOverlayImage() {
    if (!isOpen(g_orig_id)) return;
    selectImage(g_orig_id);

    roiManager("show all");

    // Flatten: burns ROI Manager overlays into a new RGB copy
    run("Flatten");
    flat_id = getImageID();

    // Set JPEG quality then save
    run("Input/Output...", "jpeg=65");
    overlay_path = out_dir + base + "_overlay.jpg";
    saveAs("Jpeg", overlay_path);
    selectImage(flat_id); close();

    // Append summary rows in canonical order:
    //   1. villus_area (sorted by VH number)
    //   2. per-loop means (VH, VW, CD, muscularis)
    //   3. VH/CD ratio
    _appendVillusAreas();
    _appendLoopAverages();
    _appendVHCDRatios();

    showMessage("Session complete",
        "Results saved to:\n  " + g_csv + "\n\n" +
        "Coordinates saved to:\n  " + g_json + "\n\n" +
        "Overlay image saved to:\n  " + overlay_path);
}

// ============================================================
//  _appendVHCDRatios
//  Reads the CSV, groups villus_height and crypt_depth values
//  by LoopPosition, then appends one VH_CD_ratio row per group.
//  Called once at Save & Exit.  Skips loop positions that lack
//  both VH and CD data.
// ============================================================
function _appendVHCDRatios() {
    if (!File.exists(g_csv)) return;
    content = File.openAsString(g_csv);
    lines   = split(content, "\n");
    n_lines = lengthOf(lines);

    // Pass 1: collect unique loop positions that have VH or CD rows.
    //         Skip VH_CD_ratio rows to avoid double-counting.
    lp_list = newArray(0);
    for (i = 1; i < n_lines; i++) {
        if (lengthOf(lines[i]) == 0) continue;
        f = split(lines[i], ",");
        if (lengthOf(f) < 6) continue;
        mtype = f[4];
        if (mtype != "villus_height" && mtype != "crypt_depth") continue;
        lp = f[3];
        found = false;
        for (j = 0; j < lengthOf(lp_list); j++) {
            if (lp_list[j] == lp) { found = true; break; }
        }
        if (!found) {
            tmp = newArray(1);  tmp[0] = lp;
            lp_list = Array.concat(lp_list, tmp);
        }
    }
    if (lengthOf(lp_list) == 0) return;

    // Pass 2: for each loop position sum VH and CD values separately.
    saved_loop = g_loop;
    for (lpi = 0; lpi < lengthOf(lp_list); lpi++) {
        lp     = lp_list[lpi];
        vh_sum = 0;  vh_n = 0;
        cd_sum = 0;  cd_n = 0;

        for (i = 1; i < n_lines; i++) {
            if (lengthOf(lines[i]) == 0) continue;
            f = split(lines[i], ",");
            if (lengthOf(f) < 6) continue;
            if (f[3] != lp) continue;
            mtype = f[4];
            val   = parseFloat(f[5]);
            if (isNaN(val)) continue;
            if (mtype == "villus_height") { vh_sum += val; vh_n++; }
            if (mtype == "crypt_depth")   { cd_sum += val; cd_n++; }
        }

        if (vh_n > 0 && cd_n > 0) {
            vh_mean = vh_sum / vh_n;
            cd_mean = cd_sum / cd_n;
            ratio   = vh_mean / cd_mean;
            ratio_notes = "nVH=" + vh_n + ";nCD=" + cd_n +
                          ";meanVH=" + d2s(vh_mean, 1) + "um" +
                          ";meanCD=" + d2s(cd_mean, 1) + "um";
            g_loop = lp;
            appendCSV("VH_CD_ratio", ratio, ratio_notes);
        }
    }
    g_loop = saved_loop;
}

// ============================================================
//  _transferToOverlayAndReset
//  Moves villus centerline ROIs (named "VH-*") from the ROI Manager
//  to the image overlay, then resets the manager.
//  Result: centerlines remain visible on-screen as reference while
//  the next measurement session starts with ROI labels from 1.
//  Polygon outlines (VW-*-outline) are discarded — too cluttered.
// ============================================================
function _transferToOverlayAndReset() {
    n = roiManager("count");
    if (n == 0) return;
    selectImage(g_orig_id);
    for (i = 0; i < n; i++) {
        roiManager("select", i);
        Overlay.addSelection();
    }
    roiManager("reset");
}

// ============================================================
//  _appendVillusAreas
//  Writes villus_area rows accumulated in g_villus_areas to the
//  CSV in VH-number order.  Called once at Save & Exit so all
//  area rows appear together after the VH/VW/CD measurement rows.
// ============================================================
function _appendVillusAreas() {
    if (lengthOf(g_villus_areas) == 0) return;
    entries = split(g_villus_areas, ";");
    saved_loop = g_loop;
    for (i = 0; i < lengthOf(entries); i++) {
        if (lengthOf(entries[i]) == 0) continue;
        parts = split(entries[i], "|");
        if (lengthOf(parts) < 3) continue;
        g_loop = parts[0];
        vn_i   = parts[1];
        area_i = parseFloat(parts[2]);
        appendCSV("villus_area", area_i, "pair=VH-" + vn_i);
    }
    g_loop = saved_loop;
}

// ============================================================
//  _appendLoopAverages
//  For each loop position that has measurements, appends one
//  mean_* row per type (villus_height, villus_width, crypt_depth,
//  muscularis).  Called once at Save & Exit.
// ============================================================
function _appendLoopAverages() {
    if (!File.exists(g_csv)) return;
    content = File.openAsString(g_csv);
    lines   = split(content, "\n");
    n_lines = lengthOf(lines);

    avg_src    = newArray("villus_height", "villus_width", "crypt_depth", "muscularis");
    avg_dst    = newArray("mean_villus_height", "mean_villus_width",
                          "mean_crypt_depth",   "mean_muscularis");

    // Collect unique loop positions
    lp_list = newArray(0);
    for (i = 1; i < n_lines; i++) {
        if (lengthOf(lines[i]) == 0) continue;
        f = split(lines[i], ",");
        if (lengthOf(f) < 6) continue;
        lp = f[3];
        found = false;
        for (j = 0; j < lengthOf(lp_list); j++) {
            if (lp_list[j] == lp) { found = true; break; }
        }
        if (!found) {
            tmp = newArray(1);  tmp[0] = lp;
            lp_list = Array.concat(lp_list, tmp);
        }
    }
    if (lengthOf(lp_list) == 0) return;

    saved_loop = g_loop;
    for (lpi = 0; lpi < lengthOf(lp_list); lpi++) {
        lp = lp_list[lpi];
        g_loop = lp;
        for (ti = 0; ti < lengthOf(avg_src); ti++) {
            src = avg_src[ti];
            sum_v = 0;  n_v = 0;
            for (i = 1; i < n_lines; i++) {
                if (lengthOf(lines[i]) == 0) continue;
                f = split(lines[i], ",");
                if (lengthOf(f) < 6) continue;
                if (f[3] != lp)  continue;
                if (f[4] != src) continue;
                val = parseFloat(f[5]);
                if (!isNaN(val)) { sum_v += val;  n_v++; }
            }
            if (n_v > 0) {
                appendCSV(avg_dst[ti], sum_v / n_v, "n=" + n_v);
            }
        }
    }
    g_loop = saved_loop;
}

// ============================================================
//  detectScaleBarPx
//  User draws a rectangle around the scale bar (black line on
//  light background).  Duplicates, converts to 8-bit, thresholds
//  dark pixels, then measures the horizontal span → pixel length
//  of the 500 um scale bar.
//  Returns the detected pixel length, or 0 on failure.
// ============================================================
function detectScaleBarPx() {
    selectImage(g_orig_id);
    setTool("rectangle");
    waitForUser("Auto-detect scale bar",
        "Draw a RECTANGLE around the scale bar line.\n" +
        "Include only the bar itself — exclude text labels.\n\n" +
        "Click OK when done.");

    if (selectionType() != 0) {
        showMessage("Cancelled", "No rectangle drawn — scale unchanged.");
        return 0;
    }

    run("Duplicate...", "title=_scale_detect");
    scale_det_id = getImageID();
    run("8-bit");

    // Threshold: select dark pixels (scale bar ≈ value 0-80 on white ≈ 200+)
    setThreshold(0, 80);
    run("Create Selection");   // creates a selection of dark pixels

    result_px = 0;
    if (selectionType() != -1) {
        getSelectionBounds(bx, by, bw, bh);
        result_px = bw;   // horizontal span = scale bar length in pixels
    } else {
        showMessage("Detection failed",
            "No dark pixels found in the selected region.\n" +
            "Try the 'Draw line' method instead.");
    }

    resetThreshold();
    run("Select None");
    selectImage(scale_det_id); close();
    selectImage(g_orig_id);
    return result_px;
}

// ============================================================
//  tuneMuscularisDetector
//  Lets the user draw a polygon around a representative muscle
//  area, then opens tune_muscularis.py in an interactive window
//  with sliders for the E and H threshold factors.
//  Saved parameters are written to musc_params.json and picked
//  up automatically by musc_thickness.py.
// ============================================================
function tuneMuscularisDetector() {
    py_tune = g_script_dir + "tune_muscularis.py";
    if (!File.exists(py_tune)) {
        showMessage("Error",
            "tune_muscularis.py not found in:\n" + g_script_dir);
        return;
    }

    selectImage(g_orig_id);
    run("Duplicate...", "title=_tune_copy");
    tune_work_id = getImageID();

    setTool("polygon");
    waitForUser("Tune: select representative muscle area",
        "Draw a POLYGON around a representative section of\n" +
        "the muscle layer (include both muscle and any crypt\n" +
        "tissue at the boundary if visible), then click OK.");

    if (selectionType() != 2) {
        if (isOpen(tune_work_id)) { selectImage(tune_work_id); close(); }
        showMessage("Cancelled", "No polygon drawn — tuning cancelled.");
        return;
    }

    run("Duplicate...", "title=_tune_crop");
    tune_crop_id = getImageID();
    getVoxelSize(tvw, tvh, tvd, tvunit);
    setBackgroundColor(255, 255, 255);
    run("Clear Outside");
    run("Select None");
    if (isOpen(tune_work_id)) { selectImage(tune_work_id); close(); }

    tmp_tune = getDirectory("temp") + "musc_tune_tmp.tif";
    selectImage(tune_crop_id);
    saveAs("Tiff", tmp_tune);
    tune_saved_id = getImageID();
    if (isOpen(tune_saved_id)) { selectImage(tune_saved_id); close(); }

    tune_cmd = g_python_cmd +
               " \"" + py_tune  + "\"" +
               " \"" + tmp_tune + "\"" +
               " "   + d2s(tvw, 8);

    // Run Python directly (no "start" wrapper).
    // "start" was tried previously but prevented matplotlib from opening its
    // GUI window.  The bare exec() call is what worked in the first test.
    // exec() may or may not block; waitForUser gives the user explicit control.
    exec("cmd.exe", "/c", tune_cmd);

    // Wait for the user to finish tuning.
    // If exec() was blocking, the slider window has already appeared and may
    // still be open (or already closed after Save & Close).
    // If exec() was non-blocking, use the taskbar to find the slider window.
    debug_log = getDirectory("temp") + "tune_muscularis_debug.txt";
    waitForUser("Muscularis Tuning — click OK when done",
        "The slider panel should be open (or may already be closed if you\n" +
        "already clicked 'Save & Close').\n\n" +
        "If you can see the panel:\n" +
        "  1. Adjust the E / H sliders until the muscle layer looks correct\n" +
        "  2. Click  'Save & Close'  in the slider panel\n\n" +
        "If the panel is not visible, check the Windows taskbar.\n\n" +
        "Click OK here to continue.\n\n" +
        "Troubleshooting — if no panel appeared at all, check:\n" +
        debug_log);

    File.delete(tmp_tune);

    // Params are now saved at [Fiji root]/GutMorphometry_config/musc_params.json
    params_file = getDirectory("imagej") + "GutMorphometry_config" +
                  File.separator + "musc_params.json";
    if (File.exists(params_file)) {
        showMessage("Tuning complete",
            "Parameters saved to:\n" + params_file + "\n\n" +
            "Future muscularis measurements will use these settings.");
    } else {
        showMessage("Tuning window closed",
            "No parameters were saved\n" +
            "(window was closed without clicking 'Save & Close').");
    }
    selectImage(g_orig_id);
}
