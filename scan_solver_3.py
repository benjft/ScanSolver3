#!/user/bin/env python3.9

"""
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- LICENCE -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Copyright © 2021 Benedict Thompson

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

                                Scan-Solver 3.0

   This program is designed to find the most efficient (time-wise) orbit for
completing a surface scan with SCAN Sat. This is done by searching orbits with
  a period p/q * T, where 'p/q' is a reduced fraction and 'T' is the sidereal
rotation period of the body being orbited, and checking for eccentricity values
               that provide complete coverage at all latitudes.

 This is done by finding values of e (between 0 and 1) for which the following
           inequality is true for all values of x (between 0 and 1):

                                   S*F >= M

                                     where

                    S = 1/x + p/q (1-e^2)^(3/2) / (1-ex)^2

is the apparent increase in fov due to latitude, and the increase in fov due to
    the planet's rotation as the satellite passes over (2pi/T / dTheta/dt)

                         F = f(a(1-e^2)/(1-ex) - R)/A

             is the field of view at each point in the orbit, and

                                  M = 180 / q

                 is the required field of view for each track.                 

                        'e' is the orbital eccentricity
                 'x' is cos of the latitude the orbit is over
           'f' is the base field of view of the scanner (in degrees)
        'a' is the semi-major axis = cube_root((p/q)^2 μ T^2 / (4 π^2))
        'μ' is the standard gravitational parameter for the body = G*m
                   'A' is the "best altitude" of the scanner

   Note that 'field of view' as used in SCAN Sat is actually the track width
 scanned, not the angle of a cone projected from the scanner. I will use it in
this way for consistency, but it is important to know when trying to understand
                          the function of each part.

       Version 1 solved circular orbits by checking that q * fov >= 180
                 Version 2 solved F >= M for elliptical orbits

       The surface component, S, added for v3.0 improves the scan times
   that can be achieved by realising that the body rotates under the orbit,
   increasing the width swept at each point. The method used keeps the maths
 simpler, but may under-estimate the gains - especially when using higher fov
                         scanners on smaller planets.
"""
from collections.abc import Iterator, Callable
from dataclasses import dataclass, field
from math import pi, sqrt, inf, gcd, ceil
from typing import Optional


@dataclass
class Body:
    """
    class for holding information about celestial bodies (planets, moons, etc)
    also includes method for calculating the semi-major axis of an orbit with
    period p/q T where T is the period of a synchronous orbit
    """
    radius: float
    rotation_period: float
    standard_gravity: float
    safe_altitude: float
    soi_radius: float

    geo_radius: float = field(init=False)

    def __post_init__(self):
        mu = self.standard_gravity
        t = self.rotation_period
        self.geo_radius = (mu * t**2 / (4 * pi**2)) ** (1/3)

    def get_sma(self, p: int, q: int) -> float:
        """
        Finds the semi-major axis required for an orbit with period (p/q)T
        where T is the sidereal rotation period of the body being orbited.
        """
        return (p/q)**(2/3) * self.geo_radius


# information on planets in the base game and outer-planets mod (prefixed opm_)
# data is from wikis so may not be 100% accurate.
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

    "eeloo": Body(210_000, 19_460, 7.4410815e10, 5_000, 119_082_940),

    # OUTER PLANETS MOD
    "opm_sarnus": Body(5_300_000, 28_500, 8.2089702e13, 580_000, 2.7401267e9),
    "opm_hale": Body(6_000, 23_555.314, 8.1199062e5, 1_000, 41_000),
    "opm_ovok": Body(26_000, 29_440.147, 1.3258591e7, 2_000, 94_000),
    "opm_eeloo": Body(210_000, 57_914.784, 7.4410815e10, 5_000, 1_159_066.2),
    "opm_slate": Body(54_000, 192_771.15, 1.9788564e12, 10_000, 9_597_157.6),
    "opm_tekto": Body(280_000, 666_154.48, 1.9244099e11, 95_000, 8_637_005.2),

    "opm_urlum": Body(2_177_000, 41_000, 1.1944574e13, 325_000, 2.5622607e9),
    "opm_polta": Body(220_000, 73_017.111, 9.0181953e10, 5_000, 1_661_114.9),
    "opm_priax": Body(74_000, 73_017.111, 3.3831766e9, 5_000, 446_767.6),
    "opm_wal": Body(370_000, 1_009_410.8, 4.9673624e11, 10_000, 18_933_505),
    "opm_tal": Body(22_000, 48_874.483, 2.1358884e8, 2_000, 139_966.65),

    "opm_neidon": Body(2_145_000, 40_250, 1.4167882e13, 260_000, 4.4157238e9),
    "opm_thatmo": Body(286_000, 306_442.67, 1.8609758e11, 35_000, 4_709_379.1),
    "opm_nissee": Body(30_000, 27_924.872, 3.9716933e8, 5_000, 7_366_476.6),

    "opm_plock": Body(189_000, 106_309.61, 5.1844895e10, 5_000, 3.1276234e8),
    "opm_karen": Body(85_050, 106_327.76, 4.6818042e9, 2_500, 939_354.32)
}


@dataclass
class Scanner:
    """Dataclass to hold scanner properties relevant to scan completion"""
    fov: float
    altitude_min: float
    altitude_best: float
    altitude_max: float


# scanners included in SCANSat (accurate as of v20.4)
SCANNERS = {
    "ms-1": Scanner(3, 20_000, 70_000, 250_000),        # Biome, VisLo
    "ms-2a": Scanner(4, 100_000, 500_000, 750_000),     # Biome, VisLo, ResLo
    "ms-r": Scanner(1.5, 70_000, 300_000, 400_000),     # Biome, VisLo, ResLo
    "r-3b": Scanner(1.5, 5_000, 70_000, 250_000),       # AltLo
    "r-eo-1": Scanner(3.5, 50_000, 100_000, 500_000),   # AltLo
    "sar-c": Scanner(3, 500_000, 700_000, 750_000),     # AltHi
    "sar-l": Scanner(4, 250_000, 500_000, 1_000_000),   # AltHi, Biome
    "sar-x": Scanner(1.5, 70_000, 250_000, 500_000),    # AltHi
    "scan-r": Scanner(1, 20_000, 70_000, 250_000),      # ResHi
    "scan-r2": Scanner(2.5, 70_000, 250_000, 500_000),  # ResHi
    "scan-rx": Scanner(3, 100_000, 500_000, 750_000),   # ResHi
    "vs-1": Scanner(1.5, 20_000, 70_000, 250_000),      # VisHi
    "vs-11": Scanner(4, 100_000, 200_000, 1_000_000),   # Anom, VisHi
    "vs-3": Scanner(2.5, 70_000, 350_000, 500_000)      # Anom, VisHi
}


@dataclass
class SolutionParams:
    """Dataclass to handle parameters for found solutions"""
    p: int
    q: int
    e_min: float
    e_max: float


# CONSTANTS USED ELSEWHERE
BOTTOM = 0
TOP = 1

FOV_MAX: float = 20  # fov capped in SCANSat to 20° after scaling
TOLERANCE: float = 1e-5
VERBOSE = False


def coprimes_of(n: int, start: int = 1, end: int = inf) -> Iterator[int]:
    """
    creates a generator to produce co-primes of n.
    :param n: The number output must be co-prime to
    :param start: The value to start at
    :param end: Maximum value to test, end after
    :return: A lazy valued iterator that outputs coprimes in range
    """
    k = start
    while k <= end:
        if gcd(n, k) == 1:
            yield k
        k += 1


def get_scaled_fov_and_altitude(scanner: Scanner, body: Body)\
        -> tuple[float, float]:
    """
    Find the maximum fov the scanner will achieve in orbit. if this exceeds
    FOV_MAX the 'ideal' altitude will be lowered until they match.
    :param scanner: The scanner being used
    :param body: The body being orbited
    :return: the fov and altitude it will be achieved at
    """
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


def find_root_near(fx: Callable[[float], float],
                   df_dx: Callable[[float], float],
                   x0: float,
                   direction: int,
                   max_dx: float = 1e-2) -> Optional[float]:
    """
    Finds a root fx(x) = 0 near x0. Only searches in the direction specified.

    Uses a modified version of the Newtonian root finding algorithm to limit
    step size and force searching in a set direction from the starting value.

    Limit to step size was required to prevent overshoot taking it too far past
    the intended root and out of the domain.


    :param fx: the function to find a root in
    :param df_dx: derivative of f(x)
    :param x0: starting value
    :param direction: direction to search in
    :param max_dx: maximum single step change in x
    :return: the 'x' value of the root
    """

    y0 = fx(x0)

    _x, x = inf, x0
    while abs(_x - x) > TOLERANCE:
        _x = x  # save old value

        y = fx(x)
        dy = df_dx(x)
        r = y/dy

        # limit distance moved
        if abs(r) > max_dx:
            r = max_dx if r > 0 else -max_dx

        # if y has flipped (+ -> - or - -> +) change direction: passed root
        sign = int(direction) if y0*y > 0 else -int(direction)

        # ensure moving in correct direction
        if r*sign < 0:
            r = -r
        x += r

        # Any x or y we pass in is limited between 0 and 1. If we exceed these
        # we risk errors. Should never happen unless there is no root.
        if not (0 <= x <= 1):
            return None
    return x


def find_root_between(fx: Callable[[float], float], x0: float, x1: float)\
        -> float:
    """
    Finds a root between the specified starting values.

    Uses a bisection algorithm to quickly find a root between the two values.
    If multiple roots exist any one of them may be found.

    Assumes that f(x0) < 0 and f(x1) > 0 so should find roots that are
    increasing when moving from x0 to x1.
    x0 does not need to be smaller than x1 for this reason. If there are no
    roots, it will more towards closest value to 0, but not guaranteed.

    :param fx: the function to find roots in
    :param x0: value where f(x) < 0
    :param x1: value were f(x) > 0
    :return: x coordinate of found root
    """

    y0, y1 = fx(x0), fx(x1)
    if VERBOSE and y0*y1 > 0:
        print("WARN: no guaranteed root between x0 and x1.")

    # bisect space until near root
    while abs(x0 - x1) > TOLERANCE:
        x = (x0 + x1) / 2
        y = fx(x)

        if y < 0:
            x0 = x
            y0 = y
        else:
            x1 = x
            y1 = y

    if y0*y1 > 0:  # both same sign -> no root found
        if VERBOSE:
            print("WARN: No root found")
        if abs(y0) < abs(y1):  # return x with smallest y (assumed closest)
            return x0
        return x1
    return (x0 + x1) / 2


def find_limit(fxy: Callable[[float, float], float],
               df_dx: Callable[[float, float], float],
               df_dy: Callable[[float, float], float],
               side: int) -> Optional[float]:
    """
    Find the root with the largest y value that can be reached from the passed
    side of the domain (BOTTOM -> y = 0, TOP -> y = 1) and return its y value.
    If there are continuous roots between, None is returned instead.
    :param fxy: The function to find roots in
    :param df_dx: dirivative in x direction
    :param df_dy: derivative in y direction
    :param side: side of plot to start on (BOTTOM: 0 or TOP: 1)
    :return: the y value of the root reached, or None if continuum
    """
    sign = 1 - 2*side  # 1 if bottom, -1 if top

    # set initial values, stats bottom right (0,1) or top left (1,0)
    x = 1 - side
    x0, x1 = x, 1-x
    y = side

    if fxy(x, y) > 0:  # no root to start from on y axis, try x axis
        y = side
        x = find_root_near(lambda _x: fxy(_x, y),
                           lambda _x: df_dx(_x, y), 1-side, -sign)
        if x is None:  # positive and no root in x -> y is valid
            return y
        x0 = x
    else:  # find starting root on y axis
        y = find_root_between(lambda _y: fxy(x, _y), side, 1 - side)

        dx = sign * df_dx(x, y)
        if dx < 0:   # implies down slope leaves domain, must be max
            return y

    # repeat until root within tolerance
    while abs(x0 - x1) > TOLERANCE:
        # decides which side of curve current root is at
        dx = sign * df_dx(x, y)
        # find root on opposing side of curve
        if dx > 0:
            x0 = x
            x1 = find_root_between(lambda _x: fxy(_x, y), x0, x1)
        else:
            x1 = x
            x0 = find_root_between(lambda _x: fxy(_x, y), x1, x0)

        # bisect found roots
        x = (x0 + x1) / 2
        z = fxy(x, y)
        # move y to new root
        y = find_root_near(lambda _y: fxy(x, _y),
                           lambda _y: df_dy(x, _y),
                           y,
                           sign if z < 0 else -sign)

        if y is None:  # no roots at this x value or error
            return None
    return y


class Solver:
    """
    Class for containing solution equations and constants. Manages a single
    scanner - body pair
    """
    def __init__(self, scanner: Scanner, body: Body):
        self.__scanner: Scanner = scanner
        self.__body: Body = body

        self.min_sma = body.radius + max(scanner.altitude_min,
                                         body.safe_altitude)

        fov, fov_alt = get_scaled_fov_and_altitude(scanner, body)
        self.fov: float = fov
        self.fov_alt: float = fov_alt
        self.k: float = 180 * self.fov_alt / self.fov

# -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= EQUATION STUFF =-=-=-=-=-=-=-=-=-=-=-=-=-=-= #

    def _s(self, p: float, q: float, x: float, y: float) -> float:
        """S component of inequality (some rearrangement done)"""
        return q * (1-x*y)**2 + p*x*(1-y*y)**(3/2)

    def _ds_dx(self, p: int, q: int, x: float, y: float) -> float:
        """partial derivative of S with respect to x"""
        return p*(1-y*y)**(3/2) - 2*q*y*(1-x*y)

    def _ds_dy(self, p: int, q: int, x: float, y: float):
        """partial derivative of S with respect to y"""
        return -3*p*x*y*(1-y**2)**(1/2) - 2*q*x*(1-y*x)

    def _f(self, p: int, q: int, x: float, y: float) -> float:
        """F component of inequality (some rearrangement done)"""
        return (1-y*y)*self.__body.get_sma(p, q) - (1-x*y)*self.__body.radius

    def _df_dx(self, p: int, q: int, x: float, y: float) -> float:
        """partial derivative of F with respect to x"""
        return self.__body.radius * y

    def _df_dy(self, p: int, q: int, x: float, y: float) -> float:
        """partial derivative of F with respect to y"""
        return x*self.__body.radius - 2*y*self.__body.get_sma(p, q)

    def _m(self, p: int, q: int, x: float, y: float) -> float:
        """M component of inequality (some rearrangement done)"""
        return self.k * x * (1 - x*y) ** 3

    def _dm_dx(self, p: int, q: int, x: float, y: float) -> float:
        """partial derivative of M with respect to x"""
        return self.k * (1 - 4*x*y) * (1 - x*y) ** 2

    def _dm_dy(self, p: int, q: int, x: float, y: float) -> float:
        """partial derivative of M with respect to y"""
        return -3 * self.k * x**2 * (1 - x*y)**2

    def _inequality_value(self, p: int, q: int, x: float, y: float) -> float:
        """Calculates difference between the two sides of the inequality"""

        s = self._s(p, q, x, y)
        f = self._f(p, q, x, y)
        m = self._m(p, q, x, y)

        return s*f - m

    def _inequality_d_dx(self, p: int, q: int, x: float, y: float) -> float:
        """gradient of inequality in x direction"""
        s = self._s(p, q, x, y)
        ds_dx = self._ds_dx(p, q, x, y)

        f = self._f(p, q, x, y)
        df_dx = self._df_dx(p, q, x, y)

        dm_dx = self._dm_dx(p, q, x, y)

        return s*df_dx + f*ds_dx - dm_dx

    def _inequality_d_dy(self, p: int, q: int, x: float, y: float) -> float:
        """gradient of inequality in y direction"""
        s = self._s(p, q, x, y)
        ds_dy = self._ds_dy(p, q, x, y)

        f = self._f(p, q, x, y)
        df_dy = self._df_dy(p, q, x, y)

        dm_dy = self._dm_dy(p, q, x, y)

        return s*df_dy + f*ds_dy - dm_dy

    def _fixed_track_value(self, p: int, q: int, x: float, y: float) -> float:
        """
        Variant of inequality with fov fixed at max. Ensures coverage above
        'best altitude' of scanner as other inequality will over scale.
        """
        s = self._s(p, q, x, y)
        m = (180/self.fov) * x * (1 - x*y)**2  # no altitude scaling
        return s - m

    def _fixed_track_d_dx(self, p: int, q: int, x: float, y: float) -> float:
        """x gradient of alternative inequality"""
        ds_dx = self._ds_dx(p, q, x, y)
        dm_dx = (180/self.fov) * (1 - x*y) * (1 - 3*x*y)  # no altitude scaling

        return ds_dx - dm_dx

    def _fixed_track_d_dy(self, p: int, q: int, x: float, y: float) -> float:
        """y gradient of alternative inequality"""
        ds_dy = self._ds_dy(p, q, x, y)
        dm_dy = -2 * (180/self.fov) * x**2 * (1 - x*y)  # no altitude scaling

        return ds_dy - dm_dy

# -=-=-=-=-=-=-=-=-=-=-=-=--=-=- SOLVER STUFF -=-=-=-=-=-=-=-=-=-=-=-=--=-=-= #

    def check_free_track(self, p: int, q: int)\
            -> Optional[tuple[float, float]]:
        """
        Checks if the original inequality with variable track width has valid y
        solutions with p & q
        :param p: the p value to test - co-prime to q
        :param q: the q value to test - co-prime to p
        :return: the found eccentricity limits, if any
        """
        a = self.__body.get_sma(p, q)
        if a < self.min_sma:
            return None

        # find the lower limit on eccentricity by increasing from bottom
        e_min = find_limit(lambda x, y: self._inequality_value(p, q, x, y),
                           lambda x, y: self._inequality_d_dx(p, q, x, y),
                           lambda x, y: self._inequality_d_dy(p, q, x, y),
                           BOTTOM)

        # none only returned if continuous < 0 from top to bottom: no solutions
        if e_min is None:
            return None

        # find the upper limit on eccentricity by descending from top
        e_max = find_limit(lambda x, y: self._inequality_value(p, q, x, y),
                           lambda x, y: self._inequality_d_dx(p, q, x, y),
                           lambda x, y: self._inequality_d_dy(p, q, x, y),
                           TOP)

        # max < min means no continuous band
        if e_max is None or e_max < e_min:
            return None

        return e_min, e_max

    def check_fixed_track(self, p: int, q: int)\
            -> Optional[tuple[float, float]]:
        """
        Solve variant equation with fixed track width. Will only be smaller
        than original above best altitude so validates solution in that
        domain.
        :param p: the p value to test - co-prime to q
        :param q: the q value to test - co-prime to p
        :return: the found eccentricity limits, if any
        """
        a = self.__body.get_sma(p, q)
        if a < self.min_sma:
            return None

        # find the lower limit on eccentricity by increasing from bottom
        e_min = find_limit(lambda x, y: self._fixed_track_value(p, q, x, y),
                           lambda x, y: self._fixed_track_d_dx(p, q, x, y),
                           lambda x, y: self._fixed_track_d_dy(p, q, x, y),
                           BOTTOM)

        # none only returned if continuous < 0 from top to bottom: no solutions
        if e_min is None:
            return None

        # find the upper limit on eccentricity by descending from top
        e_max = find_limit(lambda x, y: self._fixed_track_value(p, q, x, y),
                           lambda x, y: self._fixed_track_d_dx(p, q, x, y),
                           lambda x, y: self._fixed_track_d_dy(p, q, x, y),
                           TOP)

        # max < min means no continuous band
        if e_max is None or e_max < e_min:
            return None

        return e_min, e_max

    def get_hard_limit(self, p: int, q: int) -> Optional[tuple[float, float]]:
        """
        Checks for hard limits on eccentricity with current sma.
        These are:
         - required altitude over poles for scanning,
         - minimum altitude at periapsis to stay safe
         - maximum altitude at apoapsis to stay within scanning range and SOI
        :return: the calculated hard limits on eccentricity, or None if they
        cannot be met
        """
        body = self.__body
        scanner = self.__scanner

        a = body.get_sma(p, q)
        r = body.radius

        if a > body.soi_radius or a > (scanner.altitude_max + r):
            return None

        # altitude at periapsis
        e_safe = 1 - (r + body.safe_altitude) / a
        # altitude over poles above min altitude
        e_pole = sqrt(1 - (r + scanner.altitude_min) / a)
        # apoapsis within max altitude
        e_apo = (r + scanner.altitude_max) / a - 1
        # apoapsis within SOI
        e_soi = body.soi_radius / a - 1

        return 0, min(e_safe, e_pole, e_apo, e_soi)

    def get_validation_limits(self, p: int, q: int)\
            -> Optional[tuple[float, float]]:
        """
        get the resulting eccentricity limits from hard limits and fixed track
        """
        e_hard = self.get_hard_limit(p, q)
        if e_hard is None:
            return None

        e_limits = self.check_fixed_track(p, q)
        if e_limits is None:
            return None

        return max(e_hard[0], e_limits[0]), min(e_hard[1], e_limits[1])


def check_free_track(p: int, q: int, solvers: list[Solver])\
        -> Optional[tuple[float, float]]:
    """
    See if this is a valid solution for all solvers. return the eccentricity
    bounds
    """
    e_min, e_max = 0, 1
    # check solvers for limits
    for solver in solvers:
        e_limits = solver.check_free_track(p, q)
        if e_limits is None:  # solver has no solutions here
            return None
        e_min = max(e_limits[0], e_min)
        e_max = min(e_limits[1], e_max)
        if e_min > e_max:  # this solution not compatible with others
            return None
    return e_min, e_max


def validate_solution(solution: SolutionParams, solvers: list[Solver])\
        -> Optional[SolutionParams]:
    """
    validate the solution by applying further eccentricity limits from
    fixed track calculations or hard limits. return the updated solution, or
    None if it no longer works
    """
    e_min, e_max = solution.e_min, solution.e_max
    p, q = solution.p, solution.q

    for solver in solvers:
        e_limits = solver.get_validation_limits(p, q)
        if e_limits is None:
            return None
        e_min = max(e_min, e_limits[0])
        e_max = min(e_max, e_limits[1])

        if e_min > e_max:
            return None

    solution.e_min = e_min
    solution.e_max = e_max

    return solution


def find_fastest(body: Body, *scanners: Scanner)\
        -> Optional[list[SolutionParams]]:
    """
    Find the family of orbits that will complete the scan within the shortest
    period of time according to the equations set out in the header.

    Multiple orbits may be returned with slightly different parameters. These
    are all expected to complete in roughly the same amount of time equal to
    pT.
    (period = pT/q, q orbits required to complete scan, (pT/q) * q = pT)
    :param body: The body being orbited
    :param scanners: The scanners present on the vessel
    :return: all found orbits that complete the scan in the shortest time
    period
    """
    solvers = [Solver(scanner, body) for scanner in scanners]

    # limit on p - not sound this to be reached in practice, but better safe
    # than sorry
    limit = 360
    # if (min_fov := min(solver.fov for solver in solvers)) < 1:
    #     limit = ceil(limit / min_fov)

    # try all values of p
    for p in range(1, limit+1):
        solutions = []
        # try all co-primes of p
        for q in coprimes_of(p):
            # check for solutions with p/q
            e_limits = check_free_track(p, q, solvers)
            if e_limits is None:  # no solutions, increase p
                break
            else:  # solutions found, append and increase q
                solution = SolutionParams(p, q, e_limits[0], e_limits[1])
                solutions.append(solution)

        solutions = [validate_solution(s, solvers) for s in solutions]
        solutions = [s for s in solutions if s is not None]
        if solutions:
            return solutions
    return None


def get_user_input() -> tuple[Body, list[Scanner]]:
    CUSTOM = "custom"
    body_name = input(f"body [name|'{CUSTOM}']: ").lower()
    body = None
    if body_name == CUSTOM:
        r = float(input("\tradius [m]: "))
        t = float(input("\tperiod [s]: "))
        mu = float(input("\tmu [m3/s2]: "))
        a = float(input("\tsafe altitude [m]: "))
        soi = float(input("\tsoi radius [m]: "))
        body = Body(r, t, mu, a, soi)
    else:
        body = BODIES[body_name]

    scanner_names = input(f"scanners [name|'{CUSTOM}'](space-separated): ")\
        .lower()
    scanner_names = scanner_names.split(" ")
    scanners = []
    custom_n = 1
    for name in scanner_names:
        if name == CUSTOM:
            print(f"--custom scanner {custom_n}--")
            custom_n += 1
            f = float(input("\tfov [deg]: "))
            mn = float(input("\tmin altitude [m]: "))
            b = float(input("\tbest altitude [m]: "))
            mx = float(input("\tmax altitude [m]: "))
            scanners.append(Scanner(f, mn, b, mx))
        else:
            scanners.append(SCANNERS[name])
    return body, scanners


def main():
    body, scanners = get_user_input()

    min_alt, max_alt = 0, inf
    for scanner in scanners:
        if scanner.altitude_min + body.radius > body.soi_radius:
            print("Scanner requires altitude outside SOI")
            return
        if scanner.altitude_max < body.safe_altitude:
            print("Scanner cannot operate above minimum safe altitude")
            return
        if scanner.altitude_max < min_alt or scanner.altitude_min > max_alt:
            print("Scanner altitude bands do not overlap!")
            return
        min_alt = max(min_alt, scanner.altitude_min)
        max_alt = min(max_alt, scanner.altitude_max)

    solutions = find_fastest(body, *scanners)
    if solutions is None:
        print("No solutions found. Should not happen.")
    else:
        for solution in solutions:
            p, q = solution.p, solution.q
            e_min, e_max = solution.e_min, solution.e_max
            a = body.get_sma(p, q)

            print(f"({p:3}/{q:3})\t\ta = {a:.3f}m\t\t"
                  f"e = {e_min:.5f} to {e_max:.5f}")


if __name__ == "__main__":
    main()
