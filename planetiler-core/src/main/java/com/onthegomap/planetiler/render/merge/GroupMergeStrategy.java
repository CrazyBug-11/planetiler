package com.onthegomap.planetiler.render.merge;

import com.onthegomap.planetiler.render.TileMergeRunnable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * @author huyang
 */
public class GroupMergeStrategy implements MergeStrategy {

  @Override
  public List<TileMergeRunnable.GeometryWithTag> polygonMerge(List<TileMergeRunnable.GeometryWithTag> list,
    double totalSize, double maxSize, int z) {
    // 计算压缩比例
    double compressionRatio = Math.max(totalSize / maxSize, 1.0);
    if (compressionRatio <= 1.0) {
      return list;
    }

    // 按hash分组
    Map<String, List<TileMergeRunnable.GeometryWithTag>> groupedFeatures = list.stream()
      .collect(Collectors.groupingBy(TileMergeRunnable.GeometryWithTag::hash));

    List<TileMergeRunnable.GeometryWithTag> resultFeatures = new ArrayList<>();

    // 处理每个分组
    for (Map.Entry<String, List<TileMergeRunnable.GeometryWithTag>> group : groupedFeatures.entrySet()) {
      String groupKey = group.getKey();
      List<TileMergeRunnable.GeometryWithTag> groupFeatures = group.getValue();

      // 按面积降序排序
      groupFeatures.sort((a, b) -> Double.compare(b.area(), a.area()));

      // 计算分组总大小和阈值
      long groupTotalSize = groupFeatures.stream()
        .mapToLong(TileMergeRunnable.GeometryWithTag::size)
        .sum();
      long sizeThreshold = (long) (groupTotalSize - groupTotalSize / compressionRatio);

      // 保留大于阈值的要素
      long accumulatedSize = 0;
      int retainedCount = 0;
      for (TileMergeRunnable.GeometryWithTag feature : groupFeatures) {
        accumulatedSize += feature.size();
        if (accumulatedSize > sizeThreshold) {
          break;
        }
        retainedCount++;
      }

      // 添加未被合并的要素
      resultFeatures.addAll(groupFeatures.subList(retainedCount, groupFeatures.size()));
    }

    return resultFeatures;
  }
}
