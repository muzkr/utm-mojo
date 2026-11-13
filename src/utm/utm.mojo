import std.math as math
from std.math import sqrt, sin, cos, floor


struct Ellipsoid[dtype: DType](
    ImplicitlyCopyable
) where dtype.is_floating_point():
    """Ellipsoid parameters."""

    comptime T = Scalar[Self.dtype]

    var a: Self.T
    """Semi-major axis."""
    var b: Self.T
    """Semi-minor axis."""
    var f: Self.T
    """Flattening."""
    var e2: Self.T
    """First eccentricity squared."""
    var ep2: Self.T
    """Second eccentricity squared."""
    var n: Self.T
    """Third flattening."""

    def __init__(out self, a: Self.T, b: Self.T):
        self.a = a
        self.b = b
        self.f = (a - b) / a
        self.e2 = self.f * (2 - self.f)
        self.ep2 = self.e2 / (1 - self.e2)
        self.n = (a - b) / (a + b)

    @staticmethod
    def from_af(a: Self.T, f_inv: Self.T) -> Self:
        """Creates an ellipsoid from semi-major axis and inverse flattening.

        Args:
            a: The semi-major axis, in meters.
            f_inv: The inverse flattening `1 / f` (for example, `298.257223563`
                for WGS 84).

        Returns:
            A new `Ellipsoid` with all derived parameters computed.
        """
        return Self(a, a - a / f_inv)

    def _curv(self, phi: Self.T) -> Tuple[Self.T, Self.T]:
        """Get radius of curvature.

        Returns:
            A tuple (rho, nu):
            * rho: radius of curvature in the meridian
            * nu: radius of curvature in the prime vertical
        """

        var t1 = 1 - self.e2 * sin(phi) ** 2
        var rho = self.a * (1 - self.e2) / (t1 * sqrt(t1))
        var nu = self.a / sqrt(t1)
        return (rho, nu)

    def _S(self, phi: Self.T) -> Self.T:
        """Meridional arc."""

        var n = self.n
        var n2 = n**2
        var n3 = n2 * n
        var n4 = n2**2
        var n5 = n4 * n

        var sin_phi = sin(phi)
        var cos_phi = cos(phi)
        var sin_2phi = 2 * sin_phi * cos_phi
        var cos_2phi = cos_phi**2 - sin_phi**2
        var sin_4phi = 2 * sin_2phi * cos_2phi
        var cos_4phi = cos_2phi**2 - sin_2phi**2
        var sin_6phi = sin_4phi * cos_2phi + cos_4phi * sin_2phi
        # var cos_6phi = cos_4phi * cos_2phi - sin_4phi * sin_2phi
        var sin_8phi = 2 * sin_4phi * cos_4phi

        var A = self.a * (1 - n + (5 / 4) * (n2 - n3) + (81 / 64) * (n4 - n5))
        var B = (
            (3 / 2) * self.a * (n - n2 + (7 / 8) * (n3 - n4) + (55 / 64) * n5)
        )
        var C = (15 / 16) * self.a * (n2 - n3 + (3 / 4) * (n4 - n5))
        var D = (35 / 48) * self.a * (n3 - n4 + (11 / 16) * n5)
        var E = (315 / 512) * self.a * (n4 - n5)

        return (
            A * phi - B * sin_2phi + C * sin_4phi - D * sin_6phi + E * sin_8phi
        )


comptime wgs84[dtype: DType] where dtype.is_floating_point() = Ellipsoid[dtype](
    6_378_137.0, 6_356_752.3142
)
"""WGS84 ellipsoid parameters. Measured in meters."""


struct Utm:
    """UTM parameters."""

    comptime k0 = 0.9996
    """Central scale factor."""
    comptime FN = 10_000_000
    """False northing for the southern hemisphere. In meters."""
    comptime FE = 500_000
    """False easting. In meters."""
    comptime LAT_MIN = -80
    """Minimum valid latitude (south). In degrees."""
    comptime LAT_MAX = 84
    """Maximum valid latitude (north). In degrees."""
    comptime LON_MIN = -180
    """Minimum valid longitude (west). In degrees."""
    comptime LON_MAX = 180
    """Maximum valid longitude (east, exclusive). In degrees."""
    comptime ZONE_WIDTH = 6
    """Standard UTM zone width in degrees."""

    @staticmethod
    def _central_meridian(zone_index: Int) -> Int:
        """Central meridian in degrees."""
        return Utm.LON_MIN + Utm.ZONE_WIDTH // 2 + zone_index * Utm.ZONE_WIDTH

    @staticmethod
    def _zone_index[
        dtype: DType, //
    ](lon: Scalar[dtype]) -> Int where dtype.is_floating_point():
        return (Int(floor(lon)) - Utm.LON_MIN) / Utm.ZONE_WIDTH


def geo_to_grid[
    dtype: DType, //
](
    lat: Scalar[dtype],
    lon: Scalar[dtype],
    *,
    ellipsoid: Ellipsoid[dtype] = wgs84[dtype],
) raises -> Tuple[
    Int,
    Scalar[dtype],
    Scalar[dtype],
    Scalar[dtype],
    Scalar[dtype],
] where dtype.is_floating_point():
    """Convert geographic coordinates to UTM grid coordinates.

    Args:
        lat: Latitude in degrees. Must be in [Utm.LAT_MIN, Utm.LAT_MAX]
            (-80 to 84); beyond that UPS applies.
        lon: Longitude in degrees. Must be in [Utm.LON_MIN, Utm.LON_MAX)
            (-180 to 180).
        ellipsoid: Reference ellipsoid. Defaults to WGS84.

    Returns:
        A tuple (zone, northing, easting, convergence, scale_factor):
        * zone: UTM zone number, 1-60. Boundary meridians belong to
            the zone on their east.
        * northing: Northing in meters. Includes the false northing
            (10,000,000 m) for southern-hemisphere inputs (lat < 0).
        * easting: Easting in meters, including the 500,000 m false
            easting.
        * convergence: Meridian convergence in degrees.
        * scale_factor: Point scale factor (k0 at the central meridian).

    Raises:
        If `lat` or `lon` is out of range.

    Note:
        Standard 6-degree zones; the Norway/Svalbard zone exceptions
        are applied.
    """

    return _geo_to_grid[complete=True](lat, lon, ellipsoid=ellipsoid)


def geo_to_grid_simplified[
    dtype: DType, //
](
    lat: Scalar[dtype],
    lon: Scalar[dtype],
    *,
    ellipsoid: Ellipsoid[dtype] = wgs84[dtype],
) raises -> Tuple[
    Int,
    Scalar[dtype],
    Scalar[dtype],
] where dtype.is_floating_point():
    """Convert geographic coordinates to UTM grid coordinates.

    Like `geo_to_grid`, but omits the meridian convergence and
    point scale factor, saving their series expansions. Returns
    (zone, northing, easting); see `geo_to_grid` for argument
    and return-value semantics.
    """

    return _geo_to_grid[complete=False](lat, lon, ellipsoid=ellipsoid)


comptime _GridResult[dtype: DType, complete: Bool] = Tuple[
    Int,
    Scalar[dtype],
    Scalar[dtype],
    Scalar[dtype],
    Scalar[dtype],
] if complete else Tuple[
    Int,
    Scalar[dtype],
    Scalar[dtype],
]


def _geo_to_grid[
    dtype: DType,
    //,
    *,  # Keyword-only
    complete: Bool,
](
    lat: Scalar[dtype],
    lon: Scalar[dtype],
    *,  # Keyword-only
    ellipsoid: Ellipsoid[dtype],
) raises -> _GridResult[dtype, complete] where dtype.is_floating_point():
    if lat < Utm.LAT_MIN or lat > Utm.LAT_MAX:
        raise "invalid latitude"
    if lon < Utm.LON_MIN or lon >= Utm.LON_MAX:
        raise "invalid longitude"

    # ================================
    # UTM params

    var phi = radians(lat)
    var lamb = radians(lon)

    # 0-based!
    var zone = Utm._zone_index(lon)
    # Norway/Svalbard exceptions
    if (lat >= 56 and lat < 64) and (lon >= 3 and lon < 12):
        zone = 31
    elif lat >= 72 and lat < 84:
        if lon >= 0 and lon < 9:
            zone = 30
        elif lon >= 9 and lon < 21:
            zone = 32
        elif lon >= 21 and lon < 33:
            zone = 34
        elif lon >= 33 and lon < 42:
            zone = 36

    var lamb0 = radians(Utm._central_meridian(zone).cast[dtype]())
    var del_lamb = lamb - lamb0

    # -----------------------

    var sin_phi = sin(phi)
    var cos_phi = cos(phi)

    var cos_phi2 = cos_phi**2
    var cos_phi3 = cos_phi2 * cos_phi
    var cos_phi4 = cos_phi2**2
    var cos_phi5 = cos_phi4 * cos_phi
    var cos_phi6 = cos_phi3**2
    var cos_phi7 = cos_phi6 * cos_phi
    var cos_phi8 = cos_phi4**2

    var tan_phi = sin_phi / cos_phi
    var tan_phi2 = tan_phi**2
    var tan_phi4 = tan_phi2**2
    var tan_phi6 = tan_phi4 * tan_phi2

    var del_lamb2 = del_lamb**2
    # var del_lamb3 = del_lamb2 * del_lamb
    # var del_lamb4 = del_lamb2**2
    # var del_lamb5 = del_lamb4 * del_lamb
    # var del_lamb6 = del_lamb3**2
    # var del_lamb7 = del_lamb6 * del_lamb
    # var del_lamb8 = del_lamb4**2

    # ================================
    # Ellipsoid params

    var ep2 = ellipsoid.ep2
    var _, nu = ellipsoid._curv(phi)
    var S = ellipsoid._S(phi)

    var ep4 = ep2**2
    var ep6 = ep2 * ep4
    var ep8 = ep4**2

    # ================================

    var T1 = S * Utm.k0
    var T2 = (nu * sin_phi * cos_phi * Utm.k0) / 2
    var T3 = ((nu * sin_phi * cos_phi3 * Utm.k0) / 24) * (
        5 - tan_phi2 + 9 * ep2 * cos_phi2 + 4 * ep4 * cos_phi4
    )
    var T4 = ((nu * sin_phi * cos_phi5 * Utm.k0) / 720) * (
        61
        - 58 * tan_phi2
        + tan_phi4
        + 270 * ep2 * cos_phi2
        - 330 * tan_phi2 * ep2 * cos_phi2
        + 445 * ep4 * cos_phi4
        + 324 * ep6 * cos_phi6
        - 680 * tan_phi2 * ep4 * cos_phi4
        + 88 * ep8 * cos_phi8
        - 600 * tan_phi2 * ep6 * cos_phi6
        - 192 * tan_phi2 * ep8 * cos_phi8
    )
    var T5 = ((nu * sin_phi * cos_phi7 * Utm.k0) / 40320) * (
        1385 - 3111 * tan_phi2 + 543 * tan_phi4 - tan_phi6
    )

    var T6 = nu * cos_phi * Utm.k0
    var T7 = (nu * cos_phi3 * Utm.k0) / 6 * (1 - tan_phi2 + ep2 * cos_phi2)
    var T8 = ((nu * cos_phi5 * Utm.k0) / 120) * (
        5
        - 18 * tan_phi2
        + tan_phi4
        + 14 * ep2 * cos_phi2
        - 58 * tan_phi2 * ep2 * cos_phi2
        + 13 * ep4 * cos_phi4
        + 4 * ep6 * cos_phi6
        - 64 * tan_phi2 * ep4 * cos_phi4
        - 24 * tan_phi2 * ep6 * cos_phi6
    )
    var T9 = ((nu * cos_phi7 * Utm.k0) / 5040) * (
        61 - 479 * tan_phi2 + 179 * tan_phi4 - tan_phi6
    )

    var FN: Scalar[dtype] = Utm.FN if lat < 0 else 0

    # var N = FN + (
    #     T1 + del_lamb2 * T2 + del_lamb4 * T3 + del_lamb6 * T4 + del_lamb8 * T5
    # )
    # var E = Utm.FE + (
    #     del_lamb * T6 + del_lamb3 * T7 + del_lamb5 * T8 + del_lamb7 * T9
    # )

    # Re-write with Honor's method
    var N = FN + (
        T1
        + del_lamb2
        * (T2 + del_lamb2 * (T3 + del_lamb2 * (T4 + del_lamb2 * T5)))
    )
    var E = Utm.FE + (
        del_lamb * (T6 + del_lamb2 * (T7 + del_lamb2 * (T8 + del_lamb2 * T9)))
    )

    comptime if complete:
        var T18 = sin_phi
        var T19 = ((sin_phi * cos_phi2) / 3) * (
            1 + 3 * ep2 * cos_phi2 + 2 * ep4 * cos_phi4
        )
        var T20 = ((sin_phi * cos_phi4) / 15) * (
            2
            - tan_phi2
            + 15 * ep2 * cos_phi2
            + 35 * ep4 * cos_phi4
            - 15 * tan_phi2 * ep2 * cos_phi2
            + 33 * ep6 * cos_phi6
            - 50 * tan_phi2 * ep4 * cos_phi4
            + 11 * ep8 * cos_phi8
            - 60 * tan_phi2 * ep6 * cos_phi6
            - 24 * tan_phi2 * ep8 * cos_phi8
        )
        var T21 = ((sin_phi * cos_phi6) / 315) * (
            17 - 26 * tan_phi2 + 2 * tan_phi4
        )

        var T26 = (cos_phi2 / 2) * (1 + ep2 * cos_phi2)
        var T27 = (cos_phi4 / 24) * (
            5
            - 4 * tan_phi2
            + 14 * ep2 * cos_phi2
            + 13 * ep4 * cos_phi4
            - 28 * tan_phi2 * ep2 * cos_phi2
            + 4 * ep6 * cos_phi6
            - 48 * tan_phi2 * ep4 * cos_phi4
            - 24 * tan_phi2 * ep6 * cos_phi6
        )
        var T28 = (cos_phi6 / 720) * (61 - 148 * tan_phi2 + 16 * tan_phi4)

        # var C = (
        #     del_lamb * T18 + del_lamb3 * T19 + del_lamb5 * T20 + del_lamb7 * T21
        # )
        # var k = Utm.k0 * (
        #     1 + del_lamb2 * T26 + del_lamb4 * T27 + del_lamb6 * T28
        # )

        var C = del_lamb * (
            T18 + del_lamb2 * (T19 + del_lamb2 * (T20 + del_lamb2 * T21))
        )
        var k = Utm.k0 * (
            1 + del_lamb2 * (T26 + del_lamb2 * (T27 + del_lamb2 * T28))
        )

        return rebind[_GridResult[dtype, complete]](
            (zone + 1, N, E, degrees(C), k)
        )
    else:
        return rebind[_GridResult[dtype, complete]]((zone + 1, N, E))


def grid_to_geo[
    dtype: DType, //
](
    zone: Int,
    N: Scalar[dtype],
    E: Scalar[dtype],
    *,  # Keyword-only
    northern: Bool,
    ellipsoid: Ellipsoid[dtype] = wgs84[dtype],
) raises -> Tuple[
    Scalar[dtype],
    Scalar[dtype],
    Scalar[dtype],
    Scalar[dtype],
] where dtype.is_floating_point():
    """Convert grid coordinates to geographic coordinates.

    Args:
        zone: UTM zone number. (1-60).
        N: Northing. In meters.
        E: Easting. In meters.
        northern: Northern hemisphere?
        ellipsoid: Reference ellipsoid. Defaults to WGS84.

    Returns:
        A tuple (lat, lon, convergence, scale_factor):
        * lat: Latitude in degrees.
        * lon: Longitude in degrees.
        * convergence: Meridian convergence in degrees.
        * scale_factor: Point scale factor (k0 at the central meridian).

    Raises:
        If `zone` is out of range.
    """
    return _grid_to_geo[complete=True](
        zone, N, E, northern=northern, ellipsoid=ellipsoid
    )


def grid_to_geo_simplified[
    dtype: DType, //
](
    zone: Int,
    N: Scalar[dtype],
    E: Scalar[dtype],
    *,  # Keyword-only
    northern: Bool,
    ellipsoid: Ellipsoid[dtype] = wgs84[dtype],
) raises -> Tuple[
    Scalar[dtype],
    Scalar[dtype],
] where dtype.is_floating_point():
    """Convert grid coordinates to geographic coordinates.

    Like `grid_to_geo`, but omits the meridian convergence and
    point scale factor. Returns (lat, lon); see `grid_to_geo` for argument
    and return-value semantics.
    """
    return _grid_to_geo[complete=False](
        zone, N, E, northern=northern, ellipsoid=ellipsoid
    )


comptime _GeoResult[dtype: DType, complete: Bool] = Tuple[
    Scalar[dtype],
    Scalar[dtype],
    Scalar[dtype],
    Scalar[dtype],
] if complete else Tuple[
    Scalar[dtype],
    Scalar[dtype],
]


def _grid_to_geo[
    dtype: DType,
    //,
    *,  # Keyword-only
    complete: Bool,
](
    zone: Int,
    N: Scalar[dtype],
    E: Scalar[dtype],
    *,  # Keyword-only
    northern: Bool,
    ellipsoid: Ellipsoid[dtype],
) raises -> _GeoResult[dtype, complete] where dtype.is_floating_point():
    if zone < 1 or zone > 60:
        raise "invalid zone"

    var lamb0 = radians(Utm._central_meridian(zone - 1).cast[dtype]())

    var FN = Scalar[dtype](0 if northern else Utm.FN)

    var del_N = N - FN
    var del_E = E - Utm.FE

    # ϕ': foot point to central meridian
    # Assume a spherical Earth
    var phip = del_N / ellipsoid.a / Utm.k0
    for _ in range(8):
        var T1 = ellipsoid._S(phip) * Utm.k0

        # Tolerance: 1mm / 10mm
        comptime eps: Scalar[dtype] = 0.001 if dtype == .float64 else 0.01

        # Prevent div by zero
        if abs(T1) < eps:
            break
        if abs(T1 - del_N) < eps:
            break

        phip *= del_N / T1

    var ep2 = ellipsoid.ep2
    var rho, nu = ellipsoid._curv(phip)

    var ep4 = ep2**2
    var ep6 = ep4 * ep2
    var ep8 = ep6 * ep2

    var sin_phip = sin(phip)
    var cos_phip = cos(phip)

    var cos_phip2 = cos_phip**2
    var cos_phip4 = cos_phip2**2
    var cos_phip6 = cos_phip4 * cos_phip2
    var cos_phip8 = cos_phip6 * cos_phip2

    var tan_phip = sin_phip / cos_phip
    var tan_phip2 = tan_phip**2
    var tan_phip4 = tan_phip2**2
    var tan_phip6 = tan_phip4 * tan_phip2

    var k0: Scalar[dtype] = Utm.k0
    var k02 = k0**2
    var k03 = k02 * k0
    var k04 = k03 * k0
    var k05 = k04 * k0
    var k06 = k05 * k0
    var k07 = k06 * k0
    var k08 = k07 * k0

    # Note:
    #   del_E is normalized by nu (del_E = del_E / nu), so the T terms below carry
    #   only the k0 powers
    #

    del_E /= nu  # Note this!
    var del_E2 = del_E**2
    # var del_E3 = del_E2 * del_E
    # var del_E4 = del_E3 * del_E
    # var del_E5 = del_E4 * del_E
    # var del_E6 = del_E5 * del_E
    # var del_E7 = del_E6 * del_E
    # var del_E8 = del_E7 * del_E

    var T10 = tan_phip / (2 * rho * k02)
    var T11 = (tan_phip / (24 * rho * k04)) * (
        5
        + 3 * tan_phip2
        + ep2 * cos_phip2
        - 4 * ep4 * cos_phip4
        - 9 * tan_phip2 * ep2 * cos_phip2
    )
    var T12 = (tan_phip / (720 * rho * k06)) * (
        61
        + 90 * tan_phip2
        + 46 * ep2 * cos_phip2
        + 45 * tan_phip4
        - 252 * tan_phip2 * ep2 * cos_phip2
        - 3 * ep4 * cos_phip4
        + 100 * ep6 * cos_phip6
        - 66 * tan_phip2 * ep4 * cos_phip4
        - 90 * tan_phip4 * ep2 * cos_phip2
        + 88 * ep8 * cos_phip8
        + 225 * tan_phip4 * ep4 * cos_phip4
        + 84 * tan_phip2 * ep6 * cos_phip6
        - 192 * tan_phip2 * ep8 * cos_phip8
    )
    var T13 = (tan_phip / (40320 * rho * k08)) * (
        1385 + 3633 * tan_phip2 + 4095 * tan_phip4 + 1575 * tan_phip6
    )
    var T14 = 1 / (cos_phip * k0)
    var T15 = (1 / (6 * cos_phip * k03)) * (1 + 2 * tan_phip2 + ep2 * cos_phip2)
    var T16 = (1 / (120 * cos_phip * k05)) * (
        5
        + 6 * ep2 * cos_phip2
        + 28 * tan_phip2
        - 3 * ep4 * cos_phip4
        + 8 * tan_phip2 * ep2 * cos_phip2
        + 24 * tan_phip4
        - 4 * ep6 * cos_phip6
        + 4 * tan_phip2 * ep4 * cos_phip4
        + 24 * tan_phip2 * ep6 * cos_phip6
    )
    var T17 = (1 / (5040 * cos_phip * k07)) * (
        61 + 662 * tan_phip2 + 1320 * tan_phip4 + 720 * tan_phip6
    )

    # # Note: nu as an explicit factor
    # var phi = phip + nu * (
    #     -del_E2 * T10 + del_E4 * T11 - del_E6 * T12 + del_E8 * T13
    # )
    # var lamb = lamb0 + (
    #     del_E * T14 - del_E3 * T15 + del_E5 * T16 - del_E7 * T17
    # )

    # Note: nu as an explicit factor
    var phi = phip + nu * del_E2 * (
        -T10 + del_E2 * (T11 - del_E2 * (T12 + del_E2 * T13))
    )
    var lamb = lamb0 + del_E * (
        T14 - del_E2 * (T15 + del_E2 * (T16 - del_E2 * T17))
    )

    comptime if complete:
        var T22 = tan_phip / (k0)
        var T23 = (tan_phip / (3 * k03)) * (
            1 + tan_phip2 - ep2 * cos_phip2 - 2 * ep4 * cos_phip4
        )
        var T24 = (tan_phip / (15 * k05)) * (
            2
            + 5 * tan_phip2
            + 2 * ep2 * cos_phip2
            + 3 * tan_phip4
            + tan_phip2 * ep2 * cos_phip2
            + 9 * ep4 * cos_phip4
            + 20 * ep6 * cos_phip6
            - 7 * tan_phip2 * ep4 * cos_phip4
            - 27 * tan_phip2 * ep6 * cos_phip6
            + 11 * ep8 * cos_phip8
            - 24 * tan_phip2 * ep8 * cos_phip8
        )
        var T25 = (tan_phip / (315 * k07)) * (
            17 + 77 * tan_phip2 + 105 * tan_phip4 + 45 * tan_phip6
        )
        var T29 = (1 / (2 * k02)) * (1 + ep2 * cos_phip2)
        var T30 = (1 / (24 * k04)) * (
            1
            + 6 * ep2 * cos_phip2
            + 9 * ep4 * cos_phip4
            + 4 * ep6 * cos_phip6
            - 24 * tan_phip2 * ep4 * cos_phip4
            - 24 * tan_phip2 * ep6 * cos_phip6
        )
        var T31 = 1 / (720 * k06)

        # var C = del_E * T22 - del_E3 * T23 + del_E5 * T24 - del_E7 * T25
        # var k = k0 * (1 + del_E2 * T29 + del_E4 * T30 + del_E6 * T31)

        var C = del_E * (T22 - del_E2 * (T23 + del_E2 * (T24 - del_E2 * T25)))
        var k = k0 * (1 + del_E2 * (T29 + del_E2 * (T30 + del_E2 * T31)))

        return rebind[_GeoResult[dtype, complete]](
            (
                degrees(phi),
                degrees(lamb),
                degrees(C),
                k,
            )
        )
    else:
        return rebind[_GeoResult[dtype, complete]](
            (
                degrees(phi),
                degrees(lamb),
            )
        )


def radians[
    dtype: DType, //
](deg: Scalar[dtype]) -> Scalar[dtype] where dtype.is_floating_point():
    return (math.pi / 180) * deg


def degrees[
    dtype: DType, //
](rad: Scalar[dtype]) -> Scalar[dtype] where dtype.is_floating_point():
    return (180 / math.pi) * rad
