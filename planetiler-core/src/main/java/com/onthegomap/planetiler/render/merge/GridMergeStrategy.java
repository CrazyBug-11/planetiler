package com.onthegomap.planetiler.render.merge;

import com.onthegomap.planetiler.render.TileMergeRunnable;
import java.util.List;

/**
 * @author huyang
 */
public class GridMergeStrategy implements MergeStrategy{

  @Override
  public List<TileMergeRunnable.GeometryWithTag> polygonMerge(List<TileMergeRunnable.GeometryWithTag> list,
    double totalSize, double maxSize, int z) {
    return List.of();
  }
}
