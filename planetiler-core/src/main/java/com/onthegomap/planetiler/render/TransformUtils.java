package com.onthegomap.planetiler.render;

import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.util.AffineTransformation;

public class TransformUtils {

  public static AffineTransformation transform = new AffineTransformation();

  static {
    transform.scale(1.0 / 4096, 1.0 / 4096);
  }

  public static Geometry transform(Geometry geom) {
    return transform.transform(geom);
  }

  public static int pow2(int x) {
    int result = 1;
    for (int i = 0; i < x; i++) {
      result *= 2;
    }
    return result;
  }

  public static double division2(double x, int times) {
    for (int i = 1; i < times; i++) {
      x /= 2;
    }
    return x;
  }
}
