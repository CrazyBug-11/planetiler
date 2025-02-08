package com.onthegomap.planetiler.render.merge;

/**
 * @author huyang
 */
public class MergeFactory {

  public static MergeStrategy createMerge(String type) {
    return switch (type) {
      case "group" -> new GroupMergeStrategy();
      default -> new GridMergeStrategy();
    };
  }
}
