package com.onthegomap.planetiler.render;

import static org.locationtech.jts.geom.LinearRing.MINIMUM_VALID_SIZE;
import static org.locationtech.jts.operation.union.UnaryUnionOp.union;

import com.onthegomap.planetiler.FeatureMerge;
import com.onthegomap.planetiler.VectorTile;
import com.onthegomap.planetiler.archive.Tile;
import com.onthegomap.planetiler.archive.TileEncodingResult;
import com.onthegomap.planetiler.config.PlanetilerConfig;
import com.onthegomap.planetiler.geo.GeoUtils;
import com.onthegomap.planetiler.geo.GeometryException;
import com.onthegomap.planetiler.geo.GeometryType;
import com.onthegomap.planetiler.geo.MutableCoordinateSequence;
import com.onthegomap.planetiler.geo.TileCoord;
import com.onthegomap.planetiler.geo.VWSimplifier;
import com.onthegomap.planetiler.mbtiles.Mbtiles;
import com.onthegomap.planetiler.stats.DefaultStats;
import com.onthegomap.planetiler.util.CloseableIterator;
import com.onthegomap.planetiler.util.Gzip;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.OptionalLong;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.Collectors;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.CoordinateSequence;
import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryCollection;
import org.locationtech.jts.geom.GeometryFactory;
import org.locationtech.jts.geom.LineString;
import org.locationtech.jts.geom.LinearRing;
import org.locationtech.jts.geom.MultiLineString;
import org.locationtech.jts.geom.MultiPoint;
import org.locationtech.jts.geom.Point;
import org.locationtech.jts.geom.Polygon;
import org.locationtech.jts.geom.TopologyException;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.geom.util.GeometryFixer;
import org.locationtech.jts.index.strtree.STRtree;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class TileMergeRunnable implements Runnable {

  private static final Logger LOGGER = LoggerFactory.getLogger(TileMergeRunnable.class);

  private static final int DEFAULT_TILE_SIZE = 256;
  private static final double TILE_SCALE = 0.5d;
  public static final String LINESPACE_AREA = "Linespace_Area";
  public static final String SHAPE_AREA = "Shape_Area";
  public static final String LINESPACE_LENG = "Linespace_Leng";
  public static final String SHAPE_LENG = "Shape_Leng";

  private static final int EXTENT = 4096;

  public static final int[] gridSizeArray = new int[]{1024, 256, 128};

  private static final GridEntity[] gridEntities = new GridEntity[gridSizeArray.length];

  private static final int GRID_SIZE = 256;

  private static final double SCALE_GEOMETRY = (double) DEFAULT_TILE_SIZE / EXTENT;

  private static final int EXTENT_HALF = EXTENT / 2;

  private final TileCoord tileCoord;

  private final Mbtiles mbtiles;

  private final Mbtiles.TileWriter writer;

  private final PlanetilerConfig config;

  private final GridEntity gridEntity;

  private final double gridArea;

  Map<String, List<GeometryWithTag>> originFeatureInfos = new ConcurrentHashMap<>();

  private final GeometryFactory geometryFactory = new GeometryFactory();

  record GeometryWithTag(
    String layer,
    long id,
    Map<String, Object> tags,
    long group,
    Geometry geometry,
    int size,
    double area,
    String hash,
    double length
  ) {}

//  private static final Geometry[][] VECTOR_GRID = new Geometry[EXTENT][EXTENT];
//  private static final int GRID_WIDTH = EXTENT / GRID_SIZE;
//
//  static {
//    // 初始化一个网格
//    for (int i = 0; i < gridSizeArray.length; i++) {
//      gridEntities[i] = new GridEntity(gridSizeArray[i]);
//    }
//    List<Geometry> geometryList = new ArrayList<>();
//    GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory();
//    for (int i = 0; i < EXTENT; i += GRID_WIDTH) {
//      VECTOR_GRID[i] = new Geometry[EXTENT];
//      for (int j = 0; j < EXTENT; j += GRID_WIDTH) {
//        // 创建GeometryFactory实例
//        // 创建地理点
//        if (GRID_WIDTH == 1) {
//          Coordinate coord = new Coordinate(i, j);
//          Point point = geometryFactory.createPoint(coord);
//          VECTOR_GRID[i][j] = point;
//        } else {
//          MutableCoordinateSequence sequence = new MutableCoordinateSequence();
//          sequence.addPoint(i, j);
//          sequence.addPoint(i, j + GRID_WIDTH);
//          sequence.addPoint(i + GRID_WIDTH, j + GRID_WIDTH);
//          sequence.addPoint(i + GRID_WIDTH, j);
//          sequence.addPoint(i, j);
//          Polygon polygon = geometryFactory.createPolygon(sequence);
//          VECTOR_GRID[i][j] = polygon;
//          geometryList.add(polygon);
//        }
//      }
//    }
////    GeometryCollection gridView = geometryFactory.createGeometryCollection(
////      geometryList.toArray(new Geometry[geometryList.size()]));
////    int a = 0;
//  }

  public TileMergeRunnable(TileCoord tileCoord, Mbtiles mbtiles, Mbtiles.TileWriter writer, PlanetilerConfig config,
    GridEntity gridEntity) {
    this.tileCoord = tileCoord;
    this.mbtiles = mbtiles;
    this.writer = writer;
    this.config = config;
    this.gridEntity = gridEntity;
    this.gridArea = GeoUtils.getGridGeographicArea(tileCoord.x(), tileCoord.y(), tileCoord.z(), EXTENT,
      gridEntity.getGridWidth());
  }

  public TileCoord getTileCoord() {
    return tileCoord;
  }

  @Override
  public void run() {
    LOGGER.info("Starting tile merge thread {}", tileCoord);
    int x = tileCoord.x();
    int y = tileCoord.y();
    int z = tileCoord.z();
    int currentMaxY = (1 << z) - 1;

    // 计算子瓦片的范围 (XYZ坐标)
    int minChildX = x * 2;
    // + 1 的作用是确保计算出的范围包括了父级瓦片对应的所有子级瓦片，覆盖了父级瓦片的整个高度。 -1因为索引从零开始
    int maxChildX = x * 2 + 1;
    int minChildY = (currentMaxY - y) * 2;  // 翻转Y坐标并计算最小Y
    int maxChildY = (currentMaxY - y) * 2 + 1;  // 翻转Y坐标并计算最大Y
    try (CloseableIterator<Tile> tileIterator = mbtiles.getZoomTiles(z + 1, minChildX, minChildY, maxChildX,
      maxChildY)) {
      // 1.读取要素，拼接要素
      long totalSize = 0;
      while (tileIterator.hasNext()) {
        Tile next = tileIterator.next();
        byte[] bytes = next.bytes();
        totalSize += bytes.length;
        processTile(next, z);
      }
      VectorTile mergedTile = new VectorTile();
      VectorTile reduceTile = new VectorTile();
      long maxSize = config.maxFeatures();
      int from = 0;
      int to = 0;
      for (Map.Entry<String, List<GeometryWithTag>> entry : originFeatureInfos.entrySet()) {
        // ---------------------------V2.0---------------------------------------
        // 2. 要素像素化
        List<GeometryWithTag> geometryWithTags = gridMergeFeatures(entry.getValue(), totalSize, maxSize, z);
        // 简化要素
        geometryWithTags = simplifyGeometry(geometryWithTags);
        List<VectorTile.Feature> reduced = geometryToFeature(geometryWithTags);
        reduceTile.addLayerFeatures(entry.getKey(), reduced);
        // ---------------------------V1.0---------------------------------------
        // 2. 相同属性要素合并，不会产生空洞，不会产生数据丢失
        // 纯纯融合边界，不丢失任何数据
//        List<GeometryWithTag> merged = mergeSameFeatures(entry.getValue());
//        from += merged.size();
//        // 3. 要素重新排序、合并，目标是减少瓦片大小，根据瓦片大小进行，设置=2MB
//        List<VectorTile.Feature> reduced = mergeNearbyFeatures(merged, ratio);
//        to += reduced.size();
      }

      // 4. 一个瓦片只有一个输出结果
      byte[] encode;
      try {
        encode = reduceTile.encode();
      } catch (IllegalStateException e) {
        // TODO: 定位要素编码错误
        return;
      }
      TileEncodingResult result = new TileEncodingResult(tileCoord, Gzip.gzip(encode),
        OptionalLong.empty());
      if (encode.length > config.tileWarningSizeBytes()) {
        LOGGER.warn("{} {}kb uncompressed", tileCoord, encode.length / 1024);
      }

      LOGGER.info("[merged {}]size={} -> {}, count={} -> {}", tileCoord, totalSize, result.tileData().length, from, to);

      // 4.写入MBTiles
      // TODO: write丢失数据
      writer.write(result);
    } catch (Exception e) {
      LOGGER.error(String.format("处理瓦片 %s 时发生错误", tileCoord), e);
    }
    LOGGER.info("Ending tile merge thread {}", tileCoord);
  }

  private List<GeometryWithTag> simplifyGeometry(List<GeometryWithTag> geometryWithTags) {
    List<GeometryWithTag> simplifiedGeometries = new ArrayList<>();
    for (GeometryWithTag geometryWithTag : geometryWithTags) {
      Geometry originalGeometry = geometryWithTag.geometry();
      String geometryType = originalGeometry.getGeometryType();
      Geometry simplifiedGeometry;
      if (geometryType.equalsIgnoreCase(Geometry.TYPENAME_POLYGON) || geometryType.equalsIgnoreCase(
        Geometry.TYPENAME_MULTIPOLYGON)) {
        simplifiedGeometry = originalGeometry;

      } else if (geometryType.equalsIgnoreCase(Geometry.TYPENAME_LINESTRING) || geometryType.equalsIgnoreCase(
        Geometry.TYPENAME_MULTILINESTRING)) {
        simplifiedGeometry = new VWSimplifier().setTolerance(getDistanceTolerance()).apply(originalGeometry);
        if (simplifiedGeometry == null || !simplifiedGeometry.isValid() || simplifiedGeometry.isEmpty()) {
          if (simplifiedGeometry instanceof MultiLineString multiLineString) {
            List<LineString> result = new ArrayList<>();
            for (int i = 0; i < multiLineString.getNumGeometries(); i++) {
              LineString lineString = (LineString) multiLineString.getGeometryN(i);

              if (!lineString.isValid() || lineString.isEmpty()) {
                LineString fallbackLine = createFallbackLineString((LineString) originalGeometry.getGeometryN(i));
                if (fallbackLine != null) {
                  result.add(fallbackLine);
                }
              } else {
                result.add(lineString);
              }
            }

            if (result.isEmpty()) {
              simplifiedGeometry = handleEmptyGeometry(originalGeometry);
            } else {
              simplifiedGeometry = originalGeometry.getFactory().createMultiLineString(result.toArray(new LineString[0]));
            }

          } else if (simplifiedGeometry instanceof LineString) {
            simplifiedGeometry = handleEmptyGeometry(originalGeometry);
          }
        }

      } else {
        LOGGER.error("暂不支持 {} 类型数据简化", geometryType);
        continue;
      }

      // 创建简化后的GeometryWithTag
      simplifiedGeometries.add(
        new GeometryWithTag(
          geometryWithTag.layer(),
          geometryWithTag.id(),
          geometryWithTag.tags(),
          geometryWithTag.group(),
          simplifiedGeometry,
          geometryWithTag.size(),
          geometryWithTag.area(),
          geometryWithTag.hash(),
          geometryWithTag.length()
        )
      );
    }

    return simplifiedGeometries;
  }

  /**
   * 为空的 LineString 创建首尾两点组成的线段。
   */
  private LineString createFallbackLineString(LineString originalLineString) {
    Coordinate[] coordinates = originalLineString.getCoordinates();
    if (coordinates.length >= 2) {
      Coordinate start = coordinates[0];
      Coordinate end = coordinates[coordinates.length - 1];
      return originalLineString.getFactory().createLineString(new Coordinate[]{start, end});
    }

    return null;
  }

  /**
   * 当简化后的几何为空时，保留首尾两点组成的线段。
   */
  private Geometry handleEmptyGeometry(Geometry originalGeometry) {
    Coordinate[] coordinates = originalGeometry.getCoordinates();
    if (coordinates.length >= 2) {
      Coordinate start = coordinates[0];
      Coordinate end = coordinates[coordinates.length - 1];
      return originalGeometry.getFactory().createLineString(new Coordinate[]{start, end});
    } else {
      return originalGeometry.getFactory().createGeometryCollection(null);
    }
  }

  /**
   * 默认简化阈值
   */
  private double getDistanceTolerance() {
    int z = tileCoord.z();
    if (z >= 14) {
      return 0d;
    } else if (z >= 12) {
      return 64d;
    } else if (z >= 10) {
      return 8d;
    } else if (z >= 8) {
      return 4;
    } else if (z >= 6) {
      return 1d;
    } else if (z >= 4) {
      return 0.1d;
    } else {
      return 0d;
    }
  }


  // 拓扑验证方法保持不变
  private boolean isValidSimplifiedGeometry(Geometry original, Geometry simplified) {
    if (simplified == null || simplified.isEmpty()) {
      return false;
    }

    // 检查起点和终点
    Coordinate[] originalCoords = original.getCoordinates();
    Coordinate[] simplifiedCoords = simplified.getCoordinates();

    if (simplifiedCoords.length < 2 ||
      !originalCoords[0].equals(simplifiedCoords[0]) ||
      !originalCoords[originalCoords.length - 1].equals(simplifiedCoords[simplifiedCoords.length - 1])) {
      return false;
    }

    // 检查拓扑一致性
    Envelope originalEnvelope = original.getEnvelopeInternal();
    Envelope simplifiedEnvelope = simplified.getEnvelopeInternal();

    return simplified.isValid() &&
      (originalEnvelope.contains(simplifiedEnvelope) ||
        originalEnvelope.intersects(simplifiedEnvelope));
  }

  /**
   * 处理瓦片进行仿真变换
   *
   * @param tile
   */
  private void processTile(Tile tile, int z) throws IOException, GeometryException {
    List<VectorTile.Feature> features = VectorTile.decode(Gzip.gunzip(tile.bytes()));
    for (VectorTile.Feature feature : features) {
      try {
        // 1、拿到4096x4096要素 16MB(每个像素一个要素，会产生16MB大小)
        Geometry decode = feature.geometry().decode(1);
        // 2、将要素移动新的坐标下
        Geometry transformationGeom = simulationTransformation(decode, tile.coord());
        // 处理数据精度
//        transformationGeom = GeoUtils.snapAndFixPolygon(transformationGeom, DefaultStats.get(), "transform",
//          new PrecisionModel(1d));

        if (!transformationGeom.isEmpty()) {
          double shapeArea = getShapeArea(feature, transformationGeom);
          double shapeLength = getShapeLength(feature, transformationGeom);
          GeometryWithTag geometryWithTags =
            new GeometryWithTag(feature.layer(), feature.id(), feature.tags(),
              feature.group(),
              transformationGeom, feature.geometry().commands().length, shapeArea, feature.tags().toString(),
              shapeLength);
          originFeatureInfos.computeIfAbsent(feature.layer(), k -> new ArrayList<>())
            .add(geometryWithTags);
          LOGGER.debug("tile:{},currentZoomTiles {}", tile, originFeatureInfos.size());
        }
      } catch (AssertionError e) {
        // TODO: 定位要素解码错误
        LOGGER.error("[processTile] {}", tileCoord, e);
      }
    }
  }

  private static double getShapeArea(VectorTile.Feature feature, Geometry transformationGeom) {
    double shapeArea = transformationGeom.getArea();
    try {
      if (feature.hasTag(LINESPACE_AREA)) {
        shapeArea = Double.parseDouble(feature.getTag(LINESPACE_AREA).toString());
      }
    } catch (NumberFormatException ignore) {
    }
    try {
      if (feature.hasTag(SHAPE_AREA)) {
        shapeArea = Double.parseDouble(feature.getTag(SHAPE_AREA).toString());
      }
    } catch (NumberFormatException ignore) {
    }
    return shapeArea;
  }

  private double getShapeLength(VectorTile.Feature feature, Geometry transformationGeom) {
    double shapeLength = transformationGeom.getLength();
    try {
      if (feature.hasTag(LINESPACE_LENG)) {
        shapeLength = Double.parseDouble(feature.getTag(LINESPACE_LENG).toString());
      }
    } catch (NumberFormatException ignore) {
    }
    try {
      if (feature.hasTag(SHAPE_LENG)) {
        shapeLength = Double.parseDouble(feature.getTag(SHAPE_LENG).toString());
      }
    } catch (NumberFormatException ignore) {
    }
    return shapeLength;
  }

  /**
   * 仿真变换
   */
  private Geometry simulationTransformation(Geometry geometry, TileCoord childTile) throws GeometryException {
    int relativeX = childTile.x() & 1;
    int relativeY = childTile.y() & 1;

    AffineTransformation transform = new AffineTransformation();
    double translateX = relativeX * EXTENT * TILE_SCALE;
    double translateY = relativeY * EXTENT * TILE_SCALE;

    Geometry geom = transform.scale(TILE_SCALE, TILE_SCALE).translate(translateX, translateY).transform(geometry);
    if (geom.isEmpty()) {
      LOGGER.error("仿真变换后要素为空！");
    }

    if (!geom.isValid()) {
      geom = GeoUtils.fixPolygon(geom);
    }

    return geom;
  }

  record GeometryWithIndex(Geometry geometry, int index) {}

  private List<VectorTile.Feature> geometryToFeature(List<GeometryWithTag> list) {
    List<VectorTile.Feature> features = new ArrayList<>();
    AffineTransformation transform = new AffineTransformation();
    transform.scale(SCALE_GEOMETRY, SCALE_GEOMETRY);
    List<GeometryWithIndex> gridGeometryWithIndex = new ArrayList<>();
    for (int i = 0; i < list.size(); i++) {
      GeometryWithTag geometryWithTag = list.get(i);
      Geometry geometry = geometryWithTag.geometry();
      // 存储要素、要素ID
      gridGeometryWithIndex.add(new GeometryWithIndex(geometry, i));
    }
    // 合并同类型要素
    Map<Object, List<GeometryWithIndex>> grouped = gridGeometryWithIndex.stream()
      .collect(Collectors.groupingBy(geometryWithIndex -> list.get(geometryWithIndex.index()).tags()));
    List<GeometryWithIndex> geometries = new ArrayList<>();
    // 按照要素的标签进行合并
    grouped.forEach((o, geometryWithIndices) -> {
      if (geometryWithIndices.isEmpty()) {
        return;
      }

      Geometry geometry = GeoUtils.createGeometryCollection(
        geometryWithIndices.stream().map(GeometryWithIndex::geometry).collect(Collectors.toList()));
      Geometry union = geometry.union();
      if (union instanceof Polygon || union instanceof LineString) {
        geometries.add(new GeometryWithIndex(union, geometryWithIndices.getFirst().index()));
      } else if (union instanceof GeometryCollection collection) {
        if (union instanceof MultiLineString) {
          geometries.add(new GeometryWithIndex(union, geometryWithIndices.getFirst().index()));
        } else {
          for (int i = 0; i < collection.getNumGeometries(); i++) {
            Geometry geometryN = collection.getGeometryN(i);
            if (geometryN instanceof Point || geometryN instanceof MultiPoint) {
              // 此处如果不过滤point数据  GeometryCoordinateSequences.extractGroups会报错
              continue;
            }
            geometries.add(new GeometryWithIndex(geometryN, geometryWithIndices.getFirst().index()));
          }
        }
      }
    });

    for (GeometryWithIndex geometryWithIndex : geometries) {
      GeometryWithTag geometryWithTag = list.get(geometryWithIndex.index());
      Geometry geometry = geometryWithIndex.geometry();
      Geometry scaled = transform.transform(geometry);
      try {
        VectorTile.Feature feature = getFeature(scaled, geometryWithTag);
        if (feature != null) {
          features.add(feature);
        }
      } catch (Exception e) {
        LOGGER.error("[geometryToFeature] {}", geometryWithTag, e);
      }
    }
    return features;
  }

  private VectorTile.Feature getFeature(Geometry geometry, GeometryWithTag geometryWithTag)
    throws GeometryException {
    List<List<CoordinateSequence>> geoms = GeometryCoordinateSequences.extractGroups(geometry, 0);
    if (geoms.isEmpty()) {
      return null;
    }
    Geometry geom;
    if (geometry instanceof Polygon) {
      geom = GeometryCoordinateSequences.reassemblePolygons(geoms);
      /*
       * Use the very expensive, but necessary JTS Geometry#buffer(0) trick to repair invalid polygons (with self-
       * intersections) and JTS GeometryPrecisionReducer utility to snap polygon nodes to the vector tile grid
       * without introducing self-intersections.
       *
       * See https://docs.mapbox.com/vector-tiles/specification/#simplification for issues that can arise from naive
       * coordinate rounding.
       */
      geom = GeoUtils.snapAndFixPolygon(geom, DefaultStats.get(), "render");
      // JTS utilities "fix" the geometry to be clockwise outer/CCW inner but vector tiles flip Y coordinate,
      // so we need outer CCW/inner clockwise
      geom = geom.reverse();
    } else {
      geom = GeometryCoordinateSequences.reassembleLineStrings(geoms);
    }

    if (geom.isEmpty()) {
      return null;
    }

    return new VectorTile.Feature(geometryWithTag.layer, geometryWithTag.id,
      VectorTile.encodeGeometry(geom, 0), geometryWithTag.tags, geometryWithTag.group);
  }

  private List<GeometryWithTag> gridMergeFeatures(List<GeometryWithTag> list, long totalSize, long maxSize, int z) {
    String geometryType = list.getFirst().geometry().getGeometryType();
    if (geometryType.equalsIgnoreCase(Geometry.TYPENAME_POLYGON) || geometryType.equalsIgnoreCase(
      Geometry.TYPENAME_MULTIPOLYGON)) {
      return mergePolygon(list, (double) totalSize, (double) maxSize, z);

    } else if (geometryType.equalsIgnoreCase(Geometry.TYPENAME_LINESTRING) || geometryType.equalsIgnoreCase(
      Geometry.TYPENAME_MULTILINESTRING)) {
      return mergeLine(list);

    } else {
      // TODO 暂不支持的类型
      LOGGER.error("暂不支持 {} 类型数据像素化", geometryType);
      return Collections.emptyList();
    }
  }

  private List<GeometryWithTag> mergeLine(List<GeometryWithTag> list) {
    STRtree envelopeIndex = new STRtree();
    // 将要素添加到R树中
    for (int i = 0; i < list.size(); i++) {
      Geometry geometry = list.get(i).geometry();
      Envelope env = geometry.getEnvelopeInternal().copy();
      envelopeIndex.insert(env, i);
    }
    // 保存要素+标签
    List<GeometryWithTag> result = new ArrayList<>();

    int gridWidth = gridEntity.getGridWidth();
    Geometry[][] vectorGrid = gridEntity.getVectorGrid();
    Set<Integer> processedFeatures = new HashSet<>();

    int maxExtend = EXTENT + GridEntity.BUFFER;
    for (int i = -GridEntity.BUFFER; i < maxExtend; i += gridWidth) {
      for (int j = -GridEntity.BUFFER; j < maxExtend; j += gridWidth) {
        int arrayI = i + GridEntity.BUFFER;
        int arrayJ = j + GridEntity.BUFFER;

        // 获取网格
        Geometry geometry = vectorGrid[arrayI][arrayJ];
        if (geometry == null) {
          continue;
        }

        // 找到相交的数据
        Set<Integer> set = new HashSet<>();
        envelopeIndex.query(geometry.getEnvelopeInternal(), object -> {
          if (object instanceof Integer x) {
            // 这里使用勾股定理
            if (geometry.intersects(list.get(x).geometry())) {
              set.add(x);
            }
          }
        });
        List<Integer> sortLength = set.stream().sorted(Comparator.comparingDouble(x -> list.get(x).length())).toList();
        if (!sortLength.isEmpty()) {
          GeometryWithTag geometryWithTag = list.get(sortLength.getLast());
          // v1.0 裁剪要素会出现断裂情况
//            GeometryClipper clipper = new GeometryClipper(geometry.getEnvelopeInternal());
//            Geometry gridGeometry = clipper.clip(geometryWithTag.geometry(), false);
//
//            if (gridGeometry == null || Geometry.TYPENAME_MULTIPOINT.equalsIgnoreCase(gridGeometry.getGeometryType())
//              || Geometry.TYPENAME_POINT.equalsIgnoreCase(gridGeometry.getGeometryType()) || !gridGeometry.isValid()) {
////              gridGeometry = geometry;
//              continue;
//            }

          // v2.0
          if (processedFeatures.contains(geometryWithTag.geometry.hashCode())) {
            continue;
          }
          processedFeatures.add(geometryWithTag.geometry.hashCode());
//          Geometry gridGeometry = geometryWithTag.geometry();
//
//          GeometryWithTag resultGeom = new GeometryWithTag(geometryWithTag.layer, geometryWithTag.id,
//            geometryWithTag.tags, geometryWithTag.group, gridGeometry, geometryWithTag.size, geometryWithTag.area,
//            geometryWithTag.hash, geometryWithTag.length);
          result.add(geometryWithTag);
        }
      }
    }
    return result;
  }

  private List<GeometryWithTag> mergePolygon(List<GeometryWithTag> list, double totalSize,
    double maxSize, int z) {
    //    if (totalSize <= maxSize) {
//      return list;
//    }
    double ratio = Math.max((double) totalSize / (double) maxSize, 1.0); // 表示需要压缩多少倍数据
    int offset = Math.min((int) ratio, gridSizeArray.length - 1); // 避免超出范围
    STRtree envelopeIndex = new STRtree();
    // 将要素添加到R树中
    // TODO: 要素大小排序后，做个数量限制
    for (int i = 0; i < list.size(); i++) {
      Geometry geometry = list.get(i).geometry();
      Envelope env = geometry.getEnvelopeInternal().copy();
      envelopeIndex.insert(env, i);
    }
    // 保存要素+标签
    List<GeometryWithTag> result = new ArrayList<>();
    int gridWidth = gridEntity.getGridWidth();
    Geometry[][] vectorGrid = gridEntity.getVectorGrid();

    int maxExtend = EXTENT + GridEntity.BUFFER;
    for (int i = -GridEntity.BUFFER; i < maxExtend; i += gridWidth) {
      for (int j = -GridEntity.BUFFER; j < maxExtend; j += gridWidth) {
        int arrayI = i + GridEntity.BUFFER;
        int arrayJ = j + GridEntity.BUFFER;

        // 获取网格
        Geometry geometry = vectorGrid[arrayI][arrayJ];
        if (geometry == null) {
          continue;
        }
        // 找到相交的数据
        Set<Integer> set = new HashSet<>();
        envelopeIndex.query(geometry.getEnvelopeInternal(), object -> {
          if (object instanceof Integer x) {
            // 这里使用勾股定理
            if (geometry.isWithinDistance(list.get(x).geometry(), 0)) {
              set.add(x);
            }
          }
        });
        List<Integer> sortArea = set.stream().sorted(Comparator.comparingDouble(x -> list.get(x).area())).toList();
        if (!sortArea.isEmpty()) {
          // 找到真实面积最大的瓦片，当前像素归属到
          Integer geoIndex = sortArea.getLast();
          GeometryWithTag geometryWithTag = list.get(geoIndex);

          // 计算像素与要素的比值
          double areaRatio = gridArea / geometryWithTag.area;
          Geometry gridGeometry;
          if (areaRatio < config.rasterizeAreaThreshold() || config.rasterizeMaxZoom() - 1 == z) {
            // TODO: 处理边界裁剪 当像素膨胀率小于指定值或者当前层级为栅格化最大层级时  直接使用单位大小的像素
            gridGeometry = geometry;
            GeometryWithTag resultGeom = new GeometryWithTag(geometryWithTag.layer, geometryWithTag.id,
              geometryWithTag.tags, geometryWithTag.group, gridGeometry, geometryWithTag.size, geometryWithTag.area,
              geometryWithTag.hash, geometryWithTag.length);
            result.add(resultGeom);

          } else {
            for (Integer integer : sortArea) {
              result.add(list.get(integer));
            }
          }
        }
      }
    }
    return result;
  }

  /**
   * 纯纯融合边界，不丢失任何数据
   *
   * @param origin
   * @return
   * @throws GeometryException
   */
  private List<VectorTile.Feature> mergeSameFeatures(List<VectorTile.Feature> origin) throws GeometryException {
    if (origin.isEmpty()) {
      return origin;
    }
    if (origin.getFirst().geometry().geomType() != GeometryType.POLYGON) {
      return origin;
    }
    try {
      return FeatureMerge.mergeOverlappingPolygons(origin, 0);
    } catch (AssertionError e) {
      // TODO: 定位要素解码错误
      LOGGER.error("[processTile] {}", tileCoord, e);
      return Collections.emptyList();
    }
  }

  public record GeometryWrapper(Geometry geometry, double area, VectorTile.Feature feature, int commandLen) {}

  /**
   * 要素融合，会丢失数据
   *
   * @param origin 原始要素列表
   * @param ratio  必需大于等于1.0
   * @return
   */
  private List<VectorTile.Feature> mergeNearbyFeatures(List<VectorTile.Feature> origin, double ratio) {
    if (origin.isEmpty()) {
      return origin;
    }
    if (origin.getFirst().geometry().geomType() != GeometryType.POLYGON) {
      return origin;
    }
    int totalCommand = 0;
    for (VectorTile.Feature feature : origin) {
      totalCommand += feature.geometry().commands().length;
    }
    // TODO: 将这个参数提取到配置中，每一层都有同的参数
    double minDistAndBuffer = 1.0 / 16;
    // TODO: 将这个参数提取到配置中，每一层都有同的参数
    double minArea = 1.0 / 16;
    // 计算需要合并的要素数量
    int shrinkFeatureSize = origin.size() - (int) (origin.size() / ratio);
    int shrinkCommandSize = totalCommand - (int) (totalCommand / ratio);
    if (shrinkFeatureSize == 0) {
      // 不需要删除
      return origin;
    }
    List<GeometryWrapper> allGeometries = new ArrayList<>();
    GeometryFactory factory = new GeometryFactory();
    for (VectorTile.Feature feature : origin) {
      try {
        Geometry decode;
        try {
          decode = feature.geometry().decode();
        } catch (AssertionError e) {
          continue;
        }
        // 根据像素将矢量要素简化
        fastSimplifyFeature(feature, decode, factory, allGeometries);
      } catch (GeometryException e) {
        LOGGER.debug(e.getMessage());
      }
    }
    // TODO: 需要改成全局排序
    allGeometries.sort(Comparator.comparingDouble(GeometryWrapper::area));
    // 1、过滤出面积小于1像素（在4096*4096网格下）
    List<GeometryWrapper> small = new ArrayList<>();
    List<GeometryWrapper> middle = new ArrayList<>();
    List<GeometryWrapper> bigger = new ArrayList<>();
    int commandSize = 0;
    for (GeometryWrapper geometryWrapper : allGeometries) {
      if (commandSize < shrinkCommandSize) {
        if (geometryWrapper.area() < minArea) {
          small.add(geometryWrapper);
        } else {
          middle.add(geometryWrapper);
          commandSize += geometryWrapper.commandLen();
        }
      } else {
        bigger.add(geometryWrapper);
      }
    }
    // TODO: 将中型要素，融合到大型要素上
//    List<GeometryWrapper> merged = featureToGrid(middle);
//    bigger.addAll(merged);
    return bigger.stream().map(w -> w.feature).toList();

//    return bigger.stream()
//      .map(g -> g.feature.copyWithNewGeometry(g.geometry)).toList();
  }

  private static void fastSimplifyFeature(VectorTile.Feature feature, Geometry decode, GeometryFactory factory,
    List<GeometryWrapper> allGeometries) {
    MutableCoordinateSequence sequence = new MutableCoordinateSequence();
    Coordinate[] coordinates = decode.getCoordinates();
    if (coordinates.length == 0) {
      return;
    }
    Coordinate pre = coordinates[0];
    sequence.forceAddPoint(coordinates[0].x, coordinates[0].y);
    for (int i = 1; i < coordinates.length; i++) {
      Coordinate current = coordinates[i];
      if (!(Math.abs(current.x - pre.x) <= 1.0 / 16) || !(Math.abs(current.y - pre.y) <= 1.0 / 16)) {
        sequence.forceAddPoint(current.x, current.y);
      }
      pre = current;
    }
    sequence.forceAddPoint(coordinates[0].x, coordinates[0].y);
    if (sequence.size() < MINIMUM_VALID_SIZE) {
      return;
    }
    // 生成新的面数据
    Polygon polygon = new Polygon(new LinearRing(sequence, factory), new LinearRing[0], factory);
    allGeometries.add(
      new GeometryWrapper(polygon, Double.parseDouble(feature.getTag(SHAPE_AREA).toString()), feature,
        feature.geometry().commands().length));
  }

  private static void mergeMiddleFeatureToBigger(List<GeometryWrapper> bigger, List<GeometryWrapper> middle,
    double minDistAndBuffer) {
    // 构建一个保留要素的索引
    STRtree envelopeIndex = new STRtree();
    for (int i = 0; i < bigger.size(); i++) {
      Geometry geometry = bigger.get(i).geometry();
      Envelope env = geometry.getEnvelopeInternal().copy();
      env.expandBy(0);
      envelopeIndex.insert(env, i);
    }
    // 2、从待合并要素中遍历像素，目标是跟周围要素进行一个融合，不需要保存属性
    for (GeometryWrapper wrapper : middle) {
      Geometry geometry = wrapper.geometry();
      // 存储与待删除要素相交的要素
      Set<Integer> set = new HashSet<>();
      envelopeIndex.query(geometry.getEnvelopeInternal(), object -> {
        if (object instanceof Integer j) {
          Geometry b = bigger.get(j).geometry();
          if (geometry.isWithinDistance(b, minDistAndBuffer)) {
            set.add(j);
          }
        }
      });
      // 将需要删除的要素执行像素化分配给待合并的要素
      List<Integer> list = new ArrayList<>(set);
      list.sort(Comparator.comparingDouble(o -> bigger.get(o).area));
      if (!list.isEmpty()) {
        Integer target = list.getFirst();
        GeometryWrapper targetWrapper = bigger.get(target);
        List<Geometry> mergeList = new ArrayList<>();
        Geometry targetGeometry = targetWrapper.geometry();
        mergeList.add(targetGeometry);
        mergeList.add(geometry.buffer(minDistAndBuffer));
        // TODO: 合并并不一定会成功，需要优化合并算法
        Geometry merged = GeoUtils.createGeometryCollection(mergeList);
        try {
          merged = union(merged);
        } catch (TopologyException e) {
          // 缓冲结果有时无效，这使得 union 抛出异常，因此修复它并重试（见 #700）
          merged = GeometryFixer.fix(merged);
          merged = union(merged);
        }
        // TODO: 对合并完的要素再做一次简化算法
        targetWrapper.feature().copyWithNewGeometry(merged);
      } else {
        // TODO: 执行要素合并
      }
    }
  }


  private List<Geometry> pixelLine(Geometry line, Map<String, GeometryWithTag> gridToGeomtryMap,
    GeometryWithTag originalGeometryWithTag) {
    if (line instanceof MultiLineString multiLine) {
      for (int i = 0; i < multiLine.getNumGeometries(); i++) {
        pixelLine(multiLine.getGeometryN(i), gridToGeomtryMap, originalGeometryWithTag);
      }
      return Collections.emptyList();
    }

//    GeometryClipper clipper = new GeometryClipper(tileEnvelope);

    Coordinate[] coords = line.getCoordinates();
    if (coords.length < 2) {
      return Collections.emptyList();
    }

    Geometry[][] vectorGrid = gridEntity.getVectorGrid();
    for (int i = 0; i < coords.length - 1; i++) {
      Coordinate start = coords[i];
      Coordinate end = coords[i + 1];

      // 计算起点和终点的网格索引
      int gridWidth = gridEntity.getGridWidth();
      int startX = ((int) Math.floor(start.x + GridEntity.BUFFER)) / gridWidth * gridWidth;
      int startY = ((int) Math.floor(start.y + GridEntity.BUFFER)) / gridWidth * gridWidth;
      int endX = ((int) Math.floor(end.x + GridEntity.BUFFER)) / gridWidth * gridWidth;
      int endY = (int) Math.floor((end.y + GridEntity.BUFFER)) / gridWidth * gridWidth;

      int max = 2 * GridEntity.BUFFER;
      startX = Math.max(0, Math.min(startX, max));
      startY = Math.max(0, Math.min(startY, max));
      endX = Math.max(0, Math.min(endX, max));
      endY = Math.max(0, Math.min(endY, max));

      // 扫描起点到终点的所有网格
      for (int x = Math.min(startX, endX); x <= Math.max(startX, endX); x += gridWidth) {
        for (int y = Math.min(startY, endY); y <= Math.max(startY, endY); y += gridWidth) {
          Geometry geometry = vectorGrid[x][y];
          if (!geometry.intersects(originalGeometryWithTag.geometry)) {
            continue;
          }

          String gridKey = x + "_" + y;
          if (gridToGeomtryMap.containsKey(gridKey)) {
            GeometryWithTag existGeometry = gridToGeomtryMap.get(gridKey);
            // TODO  暂时使用要素长度计算 后续使用真实长度
            if (existGeometry.length() >= originalGeometryWithTag.length()) {
              continue;
            }
          }

          Geometry gridSquare = createGridSquare(x, y, gridWidth);
          gridToGeomtryMap.put(gridKey, new GeometryWithTag(
            originalGeometryWithTag.layer(),
            originalGeometryWithTag.id(),
            originalGeometryWithTag.tags(),
            originalGeometryWithTag.group(),
            gridSquare,
            originalGeometryWithTag.size(),
            originalGeometryWithTag.area(),
            originalGeometryWithTag.hash(),
            originalGeometryWithTag.length()
          ));
        }
      }
    }

    return Collections.emptyList();
  }

  private Geometry createGridSquare(int gridX, int gridY, int gridWidth) {
    double minX = gridX;
    double minY = gridY;
    double maxX = minX + gridWidth;
    double maxY = minY + gridWidth;

    Coordinate[] squareCoords = new Coordinate[]{
      new Coordinate(minX, minY),
      new Coordinate(maxX, minY),
      new Coordinate(maxX, maxY),
      new Coordinate(minX, maxY),
      new Coordinate(minX, minY)
    };
    return geometryFactory.createPolygon(squareCoords);
  }
}
