from dataclasses import dataclass
from math import pi, sqrt, inf, gcd, ceil


@dataclass
class Body:
    radius: float
    rotation_period: float
    standard_gravity: float
    safe_altitude: float
    soi_radius: float


BODIES = {
    "kerbol": Body(261_600_000, 432_000, 1.1723328e18, 600_000, inf),

    "moho": Body(250_000, 1_210_000, 1.6860938e11, 10_000, 9_646_663),

    "eve": Body(700_000, 80_500, 8.1717302e12, 90_000, 85_109_365),
    "gilly": Body(13_000, 28_255, 8.289_449_8e6, 5_000, 126_123.27),

    "kerbin": Body(600_000, 21_549.425, 3.5316e12, 70_000, 84_159_286),
    "mun": Body(200_000, 138_984.38, 6.5138398e10, 10_000, 2_429_559.1),
    "minmus": Body(60_000, 40_400, 1.7658e9, 10_000, 2_247_428.4),

    "duna": Body(320_000, 65_517.859, 3.0136321e11, 50_000, 47_921_949),
    "ike": Body(130_000, 65_517.862, 1.8568369e10, 10_000, 1_049_598.9),

    "dres": Body(138_000, 34_800, 2.1484489e10, 10_000, 32_832_840),

    "jool": Body(6_000_000, 36_000, 2.82528e14, 200_000, 2.4559852e9),
    "laythe": Body(500_000, 52_980.879, 1.962e12, 50_000, 3_723_645.8),
    "vall": Body(300_000, 105_962.09, 2.074815e11, 25_000, 2_406_401.4),
    "tylo": Body(600_000, 211_926.36, 2.82528e12, 30_000, 10_856_518),
    "bop": Body(65_000, 544_507.43, 2.4868349e9, 25_000, 1_221_060.9),
    "pol": Body(44_000, 901_902.62, 7.2170208e8, 5_000, 1_042_138.9),

    "eeloo": Body(210_000, 19_460, 7.4410815e10, 5_000, 1.1908294e8)
}


@dataclass
class Scanner:
    fov: float
    altitude_min: float
    altitude_best: float
    altitude_max: float


@dataclass
class SolutionParams:
    p: int
    q: int
    semi_major_axis: float
    eccentricity_min: float
    eccentricity_max: float


def coprimes_of(n: int, start: int = 1, end: int = inf):
    k = start
    while k <= end:
        if gcd(n, k) == 1:
            yield k
        k += 1


FOV_MAX = 20  # fov capped in SCANSat to 20° after scaling


def get_scaled_fov_and_altitude(scanner: Scanner, body: Body) -> tuple[float, float]:
    fov = scanner.fov
    alt = scanner.altitude_best

    r = body.radius
    r_kerbin = BODIES["kerbin"].radius

    if r < r_kerbin:  # fov ony scales for bodies smaller than kerbin
        fov *= sqrt(r_kerbin / r)

    if fov > FOV_MAX:
        alt *= FOV_MAX / fov  # lower altitude to where fov = FOV_MAX
        fov = FOV_MAX

    return fov, alt
