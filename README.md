# ScanSolver3

This is a tool created for use with Kerbal Space Program and the [`SCANsat`](https://github.com/S-C-A-N/SCANsat) mod.
It was created to answer the question

> "What orbit will complete a surface scan in the shortest time?"

to speed up mapping, and also help guarentee 100% surface coverage from the scan.

## Usage
When run, the program will prompt inputs for the planet/object being targeted, and a list of the surface scanners being used.
In the Python version any of these can be `custom` creating further prompts for the specific parameters that define it.

If multiple scanners are specified, the program finds the fastest orbit that satisfies all of them, though faster orbits may exist for individual scanners, or even each individual scanner.

The program outputs a list showing the orbital parameters required for the solutions it found (currently assumes 90 degree inclination) in the following format

```
(<surface resonance ratio>)   a = <semi-major axis value>   e = <minimum eccentricity> to <maximim eccentricity>
```

The surface resonance ratio isn't required for setting up the orbit, but can be usefull for checking everything is correct.
It has the form
```
<number of surface rotations to complete> / <number of ground tracks swept>
```
so if the incorrect number of ground tracks are shown, or the scan isn't completed after the specified number of siderial days then soemthing may be wrong with the orbit, settings, or deeper in the calculation.

The orbit should be establiushed with an argument of periapsis or 0 or 180 degrees (0 or pi radians), and the apoapsis on the sunny side of the orbit (if sunlight is required to complete any of the scans).

## Calculation
Given the equations
```
S = 1/(cos v) + (p * (1 - e^2)^(3/2) / (q * (1 - e cos(v))^2)
F = MIN[(1/A0) * (a(p, q) * (1 - e^2) / (1 - e cos(v)) - R), 1]
```
where `v` is the latiutude being scanned, `p` is the number of sidereal rotations to complete the scan in, `q` is the number of ground tracks, `e` is the orbital eccentricity, `A0` is the ideal altitude the scanner functions at (bellow this the scanner's field of view is reduced), `a` is the semi-major axis required for an orbit with surface resonance of `p/q`, and `R` is the radius of the target planet/object.\
`S` is the effective field of view increase (when looking at a cylindrical projection) due to the latitude (higher lattitudes are streached giving `1/cos v`) and the rotation of the surface under the satellite (the 2nd fraction is the surface angular speed divided by the satellite's angular speed, which when multiplied by the field of view (FoV) should give how far the planet rotated during the passage).\
`F` is the amount the FoV is decreased by due to altitude at each position in the orbit, this scales linearly between 0 and `A0`, this equiation was derrived from the equation of radius from True Anomaly. As there is no scaling beyond `A0` if this is larger than 1 we cap it.

The program finds limits on `e` for a given `p` and `q` that have solutions for all values of `v` between 0 and 90 (deg) in the inequality
```
2f * S * F > 360 (deg) / q
```
where `f` is the FoV at `A0`.\
`360/q` is the minimum FoV required to complete each ground track.

The inequality will typically look something like this

![plot of inequality](https://imgur.com/7nn3UMI.png)

where `cos v` is shown on the x axis and `e` is shown on the y axis. The program is trying to find continuous horizontal bands and return the minimum and maximum values of `e` that form them.

In the program this is done by numerially finding the roots at alternating x then y values until the maximum/minimum e valued root is found.

To use this to find the optimal solution, the program starts with `p, q = 1, 1` and increments `q` to the next value co-prime with `p` (`p/q` must be a reduced fraction) until the orbit is outside of constraints (e.g. SOI, max altitude, etc.) then increments `p` and resets `q`. The first value of `p` solutions are found for will be optimal (`p` dictates how long the scan takes to finish), and all solutions at this value of `p` are returned to allow greater flexability in setting up the orbit.
