package com.onthegomap.planetiler.render.merge;

import com.onthegomap.planetiler.render.TileMergeRunnable;
import java.util.List;

/**
 * @author huyang
 */
public interface MergeStrategy {

  /**
   * 面要素合并策略
   *
   * @param list      要素列表
   * @param totalSize 总体要素大小
   * @param maxSize   配置的大小
   * @param z         leve
   * @return
   */
  List<TileMergeRunnable.GeometryWithTag> polygonMerge(List<TileMergeRunnable.GeometryWithTag> list, double totalSize,
    double maxSize, int z);
}
