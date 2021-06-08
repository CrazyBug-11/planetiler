package com.onthegomap.flatmap.openmaptiles;

import com.onthegomap.flatmap.Arguments;
import com.onthegomap.flatmap.CommonParams;
import com.onthegomap.flatmap.FileUtils;
import com.onthegomap.flatmap.OpenMapTilesMain;
import com.onthegomap.flatmap.Translations;
import com.onthegomap.flatmap.Wikidata;
import com.onthegomap.flatmap.collections.FeatureGroup;
import com.onthegomap.flatmap.collections.FeatureSort;
import com.onthegomap.flatmap.collections.LongLongMap;
import com.onthegomap.flatmap.profiles.OpenMapTilesProfile;
import com.onthegomap.flatmap.read.NaturalEarthReader;
import com.onthegomap.flatmap.read.OpenStreetMapReader;
import com.onthegomap.flatmap.read.OsmInputFile;
import com.onthegomap.flatmap.read.ShapefileReader;
import com.onthegomap.flatmap.write.MbtilesWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class Main {

  private static final Logger LOGGER = LoggerFactory.getLogger(OpenMapTilesMain.class);

  public static void main(String[] args) throws Exception {
    Arguments arguments = Arguments.fromJvmProperties();
    var stats = arguments.getStats();
    var overallTimer = stats.startTimer("openmaptiles");
    LOGGER.info("Arguments:");
    Path sourcesDir = Path.of("data", "sources");
    OsmInputFile osmInputFile = new OsmInputFile(
      arguments.inputFile("input", "OSM input file", sourcesDir.resolve("north-america_us_massachusetts.pbf")));
    Path centerlines = arguments
      .inputFile("centerline", "lake centerlines input", sourcesDir.resolve("lake_centerline.shp.zip"));
    Path naturalEarth = arguments
      .inputFile("natural_earth", "natural earth input", sourcesDir.resolve("natural_earth_vector.sqlite.zip"));
    Path waterPolygons = arguments
      .inputFile("water_polygons", "water polygons input", sourcesDir.resolve("water-polygons-split-3857.zip"));
    Path tmpDir = arguments.file("tmpdir", "temp directory", Path.of("data", "tmp"));
    boolean fetchWikidata = arguments.get("fetch_wikidata", "fetch wikidata translations", false);
    boolean useWikidata = arguments.get("use_wikidata", "use wikidata translations", true);
    Path wikidataNamesFile = arguments.file("wikidata_cache", "wikidata cache file",
      Path.of("data", "sources", "wikidata_names.json"));
    Path mbtilesOutputPath = arguments.file("output", "mbtiles output file", Path.of("data", "massachusetts.mbtiles"));
    List<String> languages = arguments.get("name_languages", "languages to use",
      "en,ru,ar,zh,ja,ko,fr,de,fi,pl,es,be,br,he".split(","));
    CommonParams config = CommonParams.from(arguments, osmInputFile);

    if (config.forceOverwrite()) {
      FileUtils.deleteFile(mbtilesOutputPath);
    } else if (Files.exists(mbtilesOutputPath)) {
      throw new IllegalArgumentException(mbtilesOutputPath + " already exists, use force to overwrite.");
    }

    LOGGER.info("Building OpenMapTiles profile into " + mbtilesOutputPath + " in these phases:");
    if (fetchWikidata) {
      LOGGER.info("  [wikidata] Fetch OpenStreetMap element name translations from wikidata");
    }
    LOGGER.info("  [lake_centerlines] Extract lake centerlines");
    LOGGER.info("  [water_polygons] Process ocean polygons");
    LOGGER.info("  [natural_earth] Process natural earth features");
    LOGGER.info("  [osm_pass1] Pre-process OpenStreetMap input (store node locations then relation members)");
    LOGGER.info("  [osm_pass2] Process OpenStreetMap nodes, ways, then relations");
    LOGGER.info("  [sort] Sort rendered features by tile ID");
    LOGGER.info("  [mbtiles] Encode each tile and write to " + mbtilesOutputPath);

    var translations = Translations.defaultProvider(languages);
    var profile = new OpenMapTilesProfile();

    Files.createDirectories(tmpDir);
    Path nodeDbPath = tmpDir.resolve("node.db");
    LongLongMap nodeLocations = LongLongMap.newFileBackedSortedTable(nodeDbPath);
    Path featureDbPath = tmpDir.resolve("feature.db");
    FeatureSort featureDb = FeatureSort.newExternalMergeSort(tmpDir.resolve("feature.db"), config.threads(), stats);
    FeatureGroup featureMap = new FeatureGroup(featureDb, profile, stats);
    stats.monitorFile("nodes", nodeDbPath);
    stats.monitorFile("features", featureDbPath);
    stats.monitorFile("mbtiles", mbtilesOutputPath);

    if (fetchWikidata) {
      Wikidata.fetch(osmInputFile, wikidataNamesFile, config, profile, stats);
    }
    if (useWikidata) {
      translations.addTranslationProvider(Wikidata.load(wikidataNamesFile));
    }

    ShapefileReader
      .process("EPSG:3857", OpenMapTilesProfile.LAKE_CENTERLINE_SOURCE, centerlines, featureMap, config, profile,
        stats);
    ShapefileReader
      .process(OpenMapTilesProfile.WATER_POLYGON_SOURCE, waterPolygons, featureMap, config, profile, stats);
    NaturalEarthReader
      .process(OpenMapTilesProfile.NATURAL_EARTH_SOURCE, naturalEarth, tmpDir.resolve("natearth.sqlite"), featureMap,
        config, profile, stats);

    try (var osmReader = new OpenStreetMapReader(OpenMapTilesProfile.OSM_SOURCE, osmInputFile, nodeLocations, profile,
      stats)) {
      osmReader.pass1(config);
      osmReader.pass2(featureMap, config);
    }

    LOGGER.info("Deleting node.db to make room for mbtiles");
    profile.release();
    Files.delete(nodeDbPath);

    featureDb.sort();

    MbtilesWriter.writeOutput(featureMap, mbtilesOutputPath, profile, config, stats);

    overallTimer.stop();

    LOGGER.info("FINISHED!");

    stats.printSummary();

    stats.close();
  }
}