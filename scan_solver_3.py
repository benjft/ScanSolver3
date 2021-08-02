from dataclasses import dataclass


@dataclass
class Body:
    radius: float
    rotation_period: float
    standard_gravity: float
    safe_altitude: float
    soi_radius: float


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

