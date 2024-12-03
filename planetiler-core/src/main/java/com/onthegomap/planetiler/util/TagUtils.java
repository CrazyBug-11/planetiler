package com.onthegomap.planetiler.util;

import java.util.Map;
import java.util.Objects;
import java.util.Optional;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class TagUtils {
  private static final Logger logger = LoggerFactory.getLogger(TagUtils.class);

  private final static int DEFAULT_RANK = 10000;

  public static Integer getFieldValueMinZoom(Map<String, Object> tags,
    Map<String, Map<String, Integer>> zoomLevelsMap) {
    if (zoomLevelsMap == null || zoomLevelsMap.isEmpty()) {
      return null;
    }

    return zoomLevelsMap.keySet().stream()
      .map(
        zoomKey -> Optional.ofNullable(tags.get(zoomKey)).map(value -> zoomLevelsMap.get(zoomKey).get(value.toString()))
          .orElse(null))
      .filter(Objects::nonNull)
      .min(Integer::compare)
      .orElse(null);
  }

  /**
   * 值越小优先级越高
   */
  public static Integer getFieldValueRank(Map<String, Object> tags, Map<String, Map<String, Integer>> priorityLevelsMap,
    Integer rank) {
    if (priorityLevelsMap == null || priorityLevelsMap.isEmpty()) {
      return DEFAULT_RANK;
    }

    int r = rank == null ? DEFAULT_RANK : rank;
    return priorityLevelsMap.keySet().stream()
      .map(priorityKey -> Optional.ofNullable(tags.get(priorityKey))
        .map(value -> priorityLevelsMap.get(priorityKey).getOrDefault(value.toString(), r)).orElse(r))
      .min(Integer::compare)
      .orElse(r);
  }

  /**
   * 计算 rankOrder，根据 numericField 或 priorityLevelsMap，
   * numericField优先级高于priorityLevelsMap
   */
  public static Integer calculateRankOrder(Map<String, Object> tags, String sortField, boolean isAsc,
    Map<String, Map<String, Integer>> priorityLevelsMap, Integer rank) {
    if (sortField != null) {
      Object fieldValue = tags.get(sortField);
      if (fieldValue == null) {
        return DEFAULT_RANK;
      }

      int numericRank;
      try {
        numericRank = Integer.parseInt(fieldValue.toString());
        return isAsc ? numericRank : -numericRank;
      } catch (NumberFormatException e) {
        logger.error("sortField is not numeric: {}", sortField);
        return DEFAULT_RANK;
      }

    } else if (priorityLevelsMap != null) {
      // 否则使用 priorityLevelsMap 逻辑
      return getFieldValueRank(tags, priorityLevelsMap, rank);
    }
    return DEFAULT_RANK;
  }
}
