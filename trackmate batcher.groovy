import fiji.plugin.trackmate.*
import fiji.plugin.trackmate.detection.*
import fiji.plugin.trackmate.tracking.jaqaman.*
import fiji.plugin.trackmate.features.FeatureFilter
import fiji.plugin.trackmate.io.*
import fiji.plugin.trackmate.visualization.table.TrackTableView
import fiji.plugin.trackmate.gui.displaysettings.DisplaySettings

import ij.IJ
import ij.ImagePlus
import java.io.File

// --- CONFIGURATION ---
def inputDir  = "C:/Users/MQ10002204/Macquarie University/Morsch Group - Documents/Grant R/PROJECTS/HiLO_SMT_TDP43/04_MASK_MERGE"
def outputDir = "C:/Users/MQ10002204/Macquarie University/Morsch Group - Documents/Grant R/PROJECTS/HiLO_SMT_TDP43/04_MASK_MERGE"

// Detector
def CHANNEL       = 1
def RADIUS        = 0.25 as Double
def THRESHOLD     = 63.287 as Double
def SUBPIXEL      = true
def MEDIAN_FILTER = true

// Tracker
def LINK_DIST     = 0.7 as Double
def GAP_DIST      = 0.9 as Double
def MAX_FRAME_GAP = 1 as Integer

// Filters
def QUALITY_MIN   = 90 as Double
def DURATION_MIN  = 0.068  as Double

// --- BATCH LOOP ---
def outDir = new File(outputDir)
outDir.mkdirs()

def files = new File(inputDir).listFiles().findAll {
    it.name.startsWith("Merged_") && it.name.endsWith(".tif")
}

println "========================================"
println "TrackMate Batch — ${new Date().format('yyyy-MM-dd HH:mm:ss')}"
println "Input dir exists: ${new File(inputDir).exists()}"
println "Files matched: ${files.size()}"
println "========================================"

def batchStart = System.currentTimeMillis()
def completed  = 0
def failed     = 0

files.eachWithIndex { file, idx ->
    def fileStart = System.currentTimeMillis()
    println "\n[${idx+1}/${files.size()}] ${file.name}"

    def imp = IJ.openImage(file.absolutePath)
    if (imp == null) {
        println "  SKIPPED — could not open file"
        failed++
        return
    }

     // --- Build cell ROI from CH2 label mask (values 1=cyto, 2=nuc) ---
    def dupl = imp.duplicate()
    dupl.setPosition(2, 1, 1)
    def labelImp = new ImagePlus("label", dupl.getProcessor().duplicate())
    dupl.close()

    labelImp.getProcessor().setThreshold(1, 65535, ij.process.ImageProcessor.NO_LUT_UPDATE)
    IJ.run(labelImp, "Convert to Mask", "")
    IJ.run(labelImp, "Create Selection", "")
    def cellRoi = labelImp.getRoi()
    labelImp.close()

    imp.setRoi(cellRoi)

    // --- TrackMate ---
    def settings = new Settings(imp)

    settings.detectorFactory = new LogDetectorFactory()
    settings.detectorSettings = [
        'TARGET_CHANNEL'          : CHANNEL,
        'RADIUS'                  : RADIUS,
        'THRESHOLD'               : THRESHOLD,
        'DO_SUBPIXEL_LOCALIZATION': SUBPIXEL,
        'DO_MEDIAN_FILTERING'     : MEDIAN_FILTER
    ]

    settings.addSpotFilter(new FeatureFilter('QUALITY', QUALITY_MIN, true))

    settings.trackerFactory = new SimpleSparseLAPTrackerFactory()
    settings.trackerSettings = settings.trackerFactory.getDefaultSettings()
    settings.trackerSettings['LINKING_MAX_DISTANCE']     = LINK_DIST
    settings.trackerSettings['GAP_CLOSING_MAX_DISTANCE'] = GAP_DIST
    settings.trackerSettings['MAX_FRAME_GAP']            = MAX_FRAME_GAP
    settings.trackerSettings['ALLOW_GAP_CLOSING']        = true
    settings.trackerSettings['ALLOW_TRACK_SPLITTING']    = false
    settings.trackerSettings['ALLOW_TRACK_MERGING']      = false

    settings.addTrackFilter(new FeatureFilter('TRACK_DURATION', DURATION_MIN, true))
    settings.addAllAnalyzers()

    println "  Running detection + tracking..."
    def model     = new Model()
    def trackmate = new TrackMate(model, settings)

    if (!trackmate.checkInput() || !trackmate.process()) {
        println "  ERROR: ${trackmate.getErrorMessage()}"
        imp.close()
        failed++
        return
    }

    def nSpotsTotal   = model.spots.getNSpots(false)
    def nSpotsVisible = model.spots.getNSpots(true)
    def nTracks       = model.trackModel.nTracks(true)
    println "  Spots detected: ${nSpotsTotal}  |  After filters: ${nSpotsVisible}  |  Tracks: ${nTracks}"

    def baseName = file.name.replace('.tif', '')

    // --- Export allspots CSV ---
    print "  Exporting allspots CSV... "
    def csvFile = new File(outDir, baseName + "_allspots.csv")
    def fw = new java.io.FileWriter(csvFile)
    fw.write("LABEL,ID,TRACK_ID,QUALITY,POSITION_X,POSITION_Y,POSITION_Z,POSITION_T,FRAME,RADIUS,VISIBILITY,MEDIAN_INTENSITY_CH1,MEDIAN_INTENSITY_CH2\n")
    fw.write("Label,Spot ID,Track ID,Quality,X,Y,Z,T,Frame,Radius,Visibility,Median intensity ch1,Median intensity ch2,\n")
    fw.write("Label,Spot ID,Track ID,Quality,X,Y,Z,T,Frame,R,Visibility,Median ch1,Median ch2,\n")
	fw.write(",,,(quality),(microns),(microns),(microns),(sec),,(microns),,,(counts),(counts)\n")

    def spots      = model.spots
    def trackModel = model.trackModel

    def spotToTrack = [:]
    trackModel.trackIDs(true).each { tid ->
        trackModel.trackSpots(tid).each { spot ->
            spotToTrack[spot.ID()] = tid
        }
    }

     spots.iterable(false).each { spot ->
        def tid    = spotToTrack.containsKey(spot.ID()) ? spotToTrack[spot.ID()] : ""
        def label  = "ID${spot.ID()}"
        def id     = spot.ID()
        def x      = spot.getFeature('POSITION_X')
        def y      = spot.getFeature('POSITION_Y')
        def z      = spot.getFeature('POSITION_Z')
        def t      = spot.getFeature('POSITION_T')
        def frame  = spot.getFeature('FRAME').intValue()
        def qual   = spot.getFeature('QUALITY')
        def radius = spot.getFeature('RADIUS')
        def medCh1 = spot.getFeature('MEDIAN_INTENSITY_CH1')
        def medCh2 = spot.getFeature('MEDIAN_INTENSITY_CH2')
        // def medCh3 = ... ← remove this line
        fw.write("${label},${id},${tid},${qual},${x},${y},${z},${t},${frame},${radius},1,${medCh1},${medCh2}\n")
        //                                                                  removed ${medCh3} ↑
    }
    fw.close()
    println "done"

	// --- Export tracks + edges CSV (TrackMate default format) ---
    print "  Exporting tracks + edges tables... "
    def ds = new DisplaySettings()

    def tracksCsvFile = new File(outDir, baseName + "_tracks.csv")
    def edgesCsvFile  = new File(outDir, baseName + "_edges.csv")

    TrackTableView.createTrackTable(model, ds).exportToCsv(tracksCsvFile)
    TrackTableView.createEdgeTable(model, ds).exportToCsv(edgesCsvFile)
    println "done"

    // --- Export XML ---
    print "  Exporting XML... "
    def xmlFile = new File(outDir, baseName + ".xml")
    def writer  = new TmXmlWriter(xmlFile)
    writer.appendLog("Batch processed by Groovy script")
    writer.appendModel(model)
    writer.appendSettings(settings)
    writer.writeToFile()
    println "done"

    imp.close()
    completed++

    // --- Timing + ETA ---
    def fileElapsed = (System.currentTimeMillis() - fileStart) / 1000.0
    def totalElapsed = (System.currentTimeMillis() - batchStart) / 1000.0
    def avgPerFile   = totalElapsed / completed
    def remaining    = files.size() - (idx + 1)
    def etaSecs      = (avgPerFile * remaining).toInteger()
    def etaMin       = etaSecs.intdiv(60)
    def etaSec       = etaSecs % 60
    println "  Done in ${String.format('%.1f', fileElapsed)}s  |  ETA: ${etaMin}m ${etaSec}s  |  ${completed}/${files.size()} complete, ${failed} failed"
}

def totalSecs = ((System.currentTimeMillis() - batchStart) / 1000.0).toInteger()
println "\n========================================"
println "Batch complete — ${new Date().format('yyyy-MM-dd HH:mm:ss')}"
println "Completed: ${completed}  |  Failed/skipped: ${failed}"
println "Total time: ${totalSecs.intdiv(60)}m ${totalSecs % 60}s"
println "========================================"