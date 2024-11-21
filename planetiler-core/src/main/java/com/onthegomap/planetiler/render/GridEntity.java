package com.onthegomap.planetiler.render;

import com.onthegomap.planetiler.geo.MutableCoordinateSequence;
import org.geotools.geometry.jts.JTSFactoryFinder;
import org.locationtech.jts.geom.Coordinate;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.GeometryFactory;

/**
 * 全局通用网格化对象，提升性能
 */
public class GridEntity {

  private static final int EXTENT = 4096;

  public static final int BUFFER = 64;

  // 创建GeometryFactory实例
  private static final GeometryFactory geometryFactory = JTSFactoryFinder.getGeometryFactory();

  private int gridSize;

  private int gridWidth;

  private Geometry[][] vectorGrid;

  public int getGridWidth() {
    return gridWidth;
  }

  public Geometry[][] getVectorGrid() {
    return vectorGrid;
  }

  public GridEntity(int gridSize) {
    // 网格256
    this.gridSize = gridSize;
    // 一个网格宽度
    this.gridWidth = EXTENT / gridSize;
    this.vectorGrid = generateVectorGrid();
  }

  private Geometry[][] generateVectorGrid() {
    int totalSize = EXTENT + 2 * BUFFER;
    Geometry[][] grid = new Geometry[totalSize][totalSize];

    int maxExtend = EXTENT + BUFFER;
    for (int i = -BUFFER; i < maxExtend; i += gridWidth) {
      for (int j = -BUFFER; j < maxExtend; j += gridWidth) {
        int arrayI = i + BUFFER;
        int arrayJ = j + BUFFER;

        if (gridWidth == 1) {
          grid[arrayI][arrayJ] = geometryFactory.createPoint(new Coordinate(i, j));
        } else {
          // 创建多边形
          MutableCoordinateSequence sequence = new MutableCoordinateSequence();
          int left = i;
          int up = j;
          int right = Math.min(i + gridWidth, maxExtend);
          int down = Math.min(j + gridWidth, maxExtend);
          sequence.addPoint(left, up);
          sequence.addPoint(left, down);
          sequence.addPoint(right, down);
          sequence.addPoint(right, up);
          sequence.addPoint(left, up);
          if (left == right && up == down) {
            grid[arrayI][arrayJ] = geometryFactory.createPoint(new Coordinate(left, up));
          } else {
            grid[arrayI][arrayJ] = geometryFactory.createPolygon(sequence);
          }
        }
      }
    }

    return grid;
  }
}
