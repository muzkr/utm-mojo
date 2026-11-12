from std.testing import TestSuite, assert_equal, assert_almost_equal
from std.builtin.globals import global_constant
import utm


# International
comptime ellipsoid = utm.Ellipsoid[.float32](6378388.0, 6356911.94613)


def test_geo_to_grid() raises:
    _test_geo_to_grid(
        _degrees(73, 0, 0),
        _degrees(45, 0, 0),
        38,
        8_100_702.9,
        500_000.0,
        0,
        0.9996,
    )

    _test_geo_to_grid(
        _degrees(30, 0, 0),
        _degrees(102, 0, 0),  # Boundary!
        48,
        3_322_624.35,
        # 789_422.07,
        210_577.93,
        -_degrees(1, 30, 3.76),
        1.000_633_54,
    )

    _test_geo_to_grid(
        _degrees(72, 4, 32.11),
        -_degrees(113, 54, 43.321),
        12,
        8_000_000.01,
        400_000.0,
        -_degrees(2, 46, 15.31),
        0.999_722_28,
    )


def _test_geo_to_grid(
    lat: Float32,
    lon: Float32,
    zone: Int,
    N: Float32,
    E: Float32,
    C: Float32,
    k: Float32,
) raises:
    ref e = global_constant[ellipsoid]()
    var zone1, N1, E1, C1, k1 = utm.geo_to_grid(lat, lon, ellipsoid=e)

    # print(zone1, N1, E1, C1, k1)
    assert_equal(zone1, zone)
    assert_almost_equal(N1, N, atol=2)
    assert_almost_equal(E1, E, atol=2)
    assert_almost_equal(C1, C, atol=_degrees(0, 0, 1).cast[.float64]())
    assert_almost_equal(k1, k, atol=0.00001)


def test_grid_to_geo() raises:
    _test_grid_to_geo(
        48,
        3_322_824.35,
        210_577.93,
        True,
        _degrees(30, 0, 6.489),
        _degrees(101, 59, 59.805),
        -_degrees(1, 30, 4.15),
        1.000_633_54,
    )

    _test_grid_to_geo(
        47,
        3_322_824.08,
        789_411.59,
        True,
        _degrees(30, 0, 6.489),
        _degrees(101, 59, 59.805),
        _degrees(1, 30, 3.96),
        1.000_633_46,
    )

    _test_grid_to_geo(
        31,
        1_000_000,
        200_000,
        True,
        _degrees(9, 2, 10.706),
        _degrees(0, 16, 17.099),
        -_degrees(0, 25, 43.95),
        1.000_713_86,
    )

    _test_grid_to_geo(
        30,
        1_000_491.75,
        859_739.88,
        True,
        _degrees(9, 2, 10.706),
        _degrees(0, 16, 17.099),
        _degrees(0, 30, 51.72),
        1.001_201_78,
    )

    _test_grid_to_geo(
        43,
        9_000_000,
        500_000,
        True,
        _degrees(81, 3, 30.487),
        _degrees(75, 0, 0),
        0,
        0.999_6,
    )

    _test_grid_to_geo(
        30,
        4_000_000,
        700_000,
        False,
        -_degrees(54, 6, 28.992),
        _degrees(0, 3, 33.695),
        -_degrees(2, 28, 45.39),
        1.000_090_80,
    )

    _test_grid_to_geo(
        31,
        4_000_329.42,
        307_758.89,
        False,
        -_degrees(54, 6, 28.992),
        _degrees(0, 3, 33.695),
        _degrees(2, 22, 58.83),
        1.000_053_45,
    )


def _test_grid_to_geo(
    zone: Int,
    N: Float32,
    E: Float32,
    northern: Bool,
    lat: Float32,
    lon: Float32,
    C: Float32,
    k: Float32,
) raises:
    ref e = global_constant[ellipsoid]()
    var lat1, lon1, C1, k1 = utm.grid_to_geo(
        zone, N, E, northern=northern, ellipsoid=e
    )

    assert_almost_equal(lat1, lat, atol=_degrees(0, 0, 1).cast[.float64]())
    assert_almost_equal(lon1, lon, atol=_degrees(0, 0, 1).cast[.float64]())
    assert_almost_equal(C1, C, atol=_degrees(0, 0, 1).cast[.float64]())
    assert_almost_equal(k1, k, atol=0.00001)


def _degrees(d: Float32, m: Float32, s: Float32) -> Float32:
    return d + m / 60 + s / 3600


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
