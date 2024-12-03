package com.onthegomap.planetiler.render;


import static org.junit.jupiter.api.Assertions.assertEquals;

import com.onthegomap.planetiler.Planetiler;
import com.onthegomap.planetiler.config.Arguments;
import com.onthegomap.planetiler.config.PlanetilerConfig;
import com.onthegomap.planetiler.geo.TileCoord;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.ValueSource;

class TileMergeRunnableTest {

  @ParameterizedTest
  @ValueSource(strings = {
    " --zoom_levels_map={\"priority\":{\"high\":1,\"medium\":5,\"low\":9}}   "
      + " --priority_levels_map={\"priority\":{\"high\":0,\"medium\":4,\"low\":9}}  "
  })
  public void testSortListByZoomAndRankAndLength(String args) {
    List<TileMergeRunnable.GeometryWithTag> list = List.of(
      new TileMergeRunnable.GeometryWithTag("layer1", 1, Map.of("priority", "high"), 1L, null, 0, 0, "hash1", 10.0),
      new TileMergeRunnable.GeometryWithTag("layer2", 2, Map.of("priority", "medium"), 1L, null, 0, 0, "hash2", 8.0),
      new TileMergeRunnable.GeometryWithTag("layer3", 3, Map.of("priority", "high"), 1L, null, 0, 0, "hash3", 15.0),
      new TileMergeRunnable.GeometryWithTag("layer4", 4, Map.of("priority", "low"), 1L, null, 0, 0, "hash4", 5.0)
    );

    Set<Integer> set = Set.of(0, 1, 2, 3);
    PlanetilerConfig config = Planetiler.create(Arguments.fromArgs((args).split("\\s+"))).config();
    TileMergeRunnable runnable = new TileMergeRunnable(TileCoord.ofXYZ(1, 1, 6), null, null, config, null);
    List<Integer> sortedList = runnable.getSortList(list, set);
    assertEquals(List.of(1, 0, 2), sortedList);

    TileMergeRunnable runnable1 = new TileMergeRunnable(TileCoord.ofXYZ(1, 1, 9), null, null, config, null);
    List<Integer> sortedList1 = runnable1.getSortList(list, set);
    assertEquals(List.of(3, 1, 0, 2), sortedList1);
  }

  @ParameterizedTest
  @ValueSource(strings = {
    " --priority_levels_map={\"priority\":{\"high\":0,\"medium\":4,\"low\":9}}  "
      + " --sort_field=num  --is_asc=true "
  })
  public void testSortDescendingOrder(String args) {
    List<TileMergeRunnable.GeometryWithTag> list = List.of(
      new TileMergeRunnable.GeometryWithTag("layer1", 1, Map.of("priority", "high", "num", "9"),
        1L, null, 0, 0, "hash1", 10.0),
      new TileMergeRunnable.GeometryWithTag("layer2", 2, Map.of("priority", "medium", "num", "1"),
        1L, null, 0, 0, "hash2", 8.0),
      new TileMergeRunnable.GeometryWithTag("layer3", 3, Map.of("priority", "low", "num", "0"),
        1L, null, 0, 0, "hash3", 15.0),
      new TileMergeRunnable.GeometryWithTag("layer4", 4, Map.of("priority", "low"),
        1L, null, 0, 0, "hash4", 15.0),
      new TileMergeRunnable.GeometryWithTag("layer5", 5, Map.of("priority", "medium", "num", "1"),
        1L, null, 0, 0, "hash5", 15.0)
    );

    Set<Integer> set = Set.of(0, 1, 2, 3, 4);

    PlanetilerConfig config = Planetiler.create(Arguments.fromArgs((args).split("\\s+"))).config();
    TileMergeRunnable runnable = new TileMergeRunnable(TileCoord.ofXYZ(1, 1, 6), null, null, config, null);
    List<Integer> sortedList = runnable.getSortList(list, set);
    assertEquals(List.of(2, 4, 1, 0, 3), sortedList);
  }
}
