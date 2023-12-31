// ################################# LICENCE ##################################
//                      Copyright © 2021 Benedict Thompson                     
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
//   of this software and associated documentation files (the "Software"), to  
// deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
//  sell copies of the Software, and to permit persons to whom the Software is 
//           furnished to do so, subject to the following conditions:          
// 
//  The above copyright notice and this permission notice shall be included in 
//             all copies or substantial portions of the Software.            
// 
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
//   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,  
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER   
//   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING  
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
//                               IN THE SOFTWARE.                              
// ############################################################################
// 
//                               Scan-Solver 3.0                              
// 
//  This program is designed to find the most efficient (time-wise) orbit for 
//  completing a surface scan with SCAN Sat. This is done by searching orbits 
//   with a period rp/rq * T, where 'rp/rq' is a reduced fraction and 'T' is the  
//     sidereal rotation period of the body being orbited, and checking for    
//     eccentricity values that provide complete coverage at all latitudes.    
// 
//     This is done by finding values of e (between 0 and 1) for which the    
//     following inequality is true for all values of x (between 0 and 1):    
// 
//                                   S*F >= M                                  
// 
//                                    where                                   
// 
//                    S = 1/x + rp/rq (1-e^2)^(3/2) / (1-ex)^2                   
// 
// is the apparent increase in fov due to latitude, and the increase in fov due
//  to the planet's rotation as the satellite passes over (2pi/T / dTheta/dt) 
// 
//                         F = f(a(1-e^2)/(1-ex) - rr)/A                        
// 
//             is the field of view at each point in the orbit, and            
// 
//                                 M = 180 / rq                                
// 
//                is the required field of view for each track.                
// 
//                       'e' is the orbital eccentricity                      
//                 'x' is cos of the latitude the orbit is over                
//          'f' is the base field of view of the scanner (in degrees)         
//       'a' is the semi-major axis = cube_root((rp/rq)^2 μ T^2 / (4 π^2))      
//        'μ' is the standard gravitational parameter for the body = G*m       
//                  'A' is the "best altitude" of the scanner                 
// 
//  Note that 'field of view' as used in SCAN Sat is actually the track width 
//  scanned, not the angle of a cone projected from the scanner. I will use it 
//   in this way for consistency, but it is important to know when trying to  
//                    understand the function of each part.                   
// 
//       Version 1 solved circular orbits by checking that rq * fov >= 180      
//                Version 2 solved F >= M for elliptical orbits               
// 
// The surface component, S, added for v3.0 improves the scan times that can be
// achieved by realising that the body rotates under the orbit, increasing the
// width swept at each point. The method used keeps the maths simpler, but may
//   under-estimate the gains - especially when using higher fov scanners on  
//                               smaller planets.                              


@lazyGlobal off.

// ################################ CONSTANTS #################################
local PI is constant:pi.
local FOV_MAX is 20.
local MAX_STEP is 1e-2.
local TOLERANCE is 1e-5.
local NONE is "NONE".

local DEBUG_ROW is -1.
// ################################ STRUCTURES ################################
// Create a new scanner structure. Dataclass for organisation.
function NewScanner {
    parameter fov, alt_min, alt_best, alt_max.
    return lex(
        "fov", fov,
        "alt_min", alt_min,
        "alt_best", alt_best,
        "alt_max", alt_max
    ).
}

// Create a new solutions params structure. Dataclass for organisation.
function NewSolutionParams {
    parameter rp, rq, sma, e_min, e_max.

    return lex(
        "rp", rp, 
        "rq", rq, 
        "sma", sma, 
        "e_min", e_min, 
        "e_max", e_max
    ).
}

// Create a new equation holder structure.
function NewEquationHolder {
    parameter target_body, alt_safe, scanner.

    local sma_min is target_body:radius + max(alt_safe, scanner:alt_min).
    local sma_max is min(target_body:radius + scanner:alt_max, 
                         target_body:soiRadius).
    local fov_alt is scaleFov(target_body, scanner).
    
    local fov is fov_alt[0].
    local alt_best is fov_alt[1].

    local k is 180 * alt_best / fov.
    local k_fixed is 180 / fov.
    
    local sma_sync is (target_body:mu * target_body:rotationperiod^2 / (4 * PI^2)) ^ (1/3).

    return lex(
        "body", target_body,
        "body_radius", target_body:radius,
        "scanner", scanner,
        "sma_min", sma_min,
        "sma_max", sma_max,
        "alt_safe", alt_safe,
        "fov", fov,
        "alt_best", alt_best,
        "k", k,
        "k_fixed", k_fixed,
        "sma_sync", sma_sync,
        "getSma", {parameter rp, rq. return sma_sync * (rp/rq)^(2/3).}
    ).
}

// ######################### Equation Holder Methods ##########################

// S component of inequality (some rearrangement done)
function EH_Surface {
    parameter rp, rq, x, y.
    return rq*(1 - x*y)^2 + rp*x*(1 - y^2)^(1.5).
}

// partial derivative of S with respect to x
function EH_SurfaceDx {
    parameter rp, rq, x, y.
    return rp*(1 - y^2)^(1.5) - 2*rq*y*(1 - x*y).
}

// partial derivative of S with respect to y
function EH_SurfaceDy {
    parameter rp, rq, x, y.
    return -3*rp*x*y*sqrt(1 - y^2) - 2*rq*x*(1 - y*x).
}

// F component of inequality (some rearrangement done)
function EH_Fov {
    parameter self, rp, rq, x, y.
    return (1 - y^2)*self:getSma(rp, rq) - (1 - x*y)*self:body_radius.
}

// partial derivative of F with respect to x
function EH_FovDx {
    parameter self, y.
    return self:body_radius*y.
}

// partial derivative of F with respect to y
function EH_FovDy {
    parameter self, rp, rq, x, y.
    return x*self:body_radius - 2*y*self:getSma(rp, rq).
}

// M component of inequality (some rearrangement done)
function EH_Min {
    parameter self, x, y.
    return self:k*x*(1 - x*y)^3.
}

// partial derivative of M with respect to x
function EH_MinDx {
    parameter self, x, y.
    return self:k*(1 - 4*x*y) * (1 - x*y)^2.
}

// partial derivative of M with respect to y
function EH_MinDy {
    parameter self, x, y.
    return -3 * self:k * x^2 * (1 - x*y)^2.
}

// Calculates difference between the two sides of the inequality
function EH_FreeTrack {
    parameter self, rp, rq, x, y.

    local s is EH_Surface(rp, rq, x, y).
    local f is EH_Fov(self, rp, rq, x, y).
    local m is EH_Min(self, x, y).

    return s * f - m.
}

// gradient of inequality in x direction
function EH_FreeTrackDx {
    parameter self, rp, rq, x, y.
    local s is EH_Surface(rp, rq, x, y).
    local ds_dx is EH_SurfaceDx(rp, rq, x, y).

    local f is EH_Fov(self, rp, rq, x, y).
    local df_dx is EH_FovDx(self, y).

    local dm_dx is EH_MinDx(self, x, y).

    return s * df_dx + f * ds_dx - dm_dx.
}

// gradient of inequality in y direction
function EH_FreeTrackDy {
    parameter self, rp, rq, x, y.
    local s is EH_Surface(rp, rq, x, y).
    local ds_dy is EH_SurfaceDy(rp, rq, x, y).

    local f is EH_Fov(self, rp, rq, x, y).
    local df_dy is EH_FovDy(self, rp, rq, x, y).

    local dm_dy is EH_MinDy(self, x, y).

    return s * df_dy + f * ds_dy - dm_dy.
}

// Variant of inequality with fov fixed at max. Ensures coverage above
// 'best altitude' of scanner as other inequality will over scale.
function EH_FixedTrack {
    parameter self, rp, rq, x, y.
    local s is EH_Surface(rp, rq, x, y).
    local m is self:k_fixed * x * (1 - x*y)^2.  // no altitude scaling
    return s - m.
}

// x gradient of alternative inequality
function EH_FixedTrackDx {
    parameter self, rp, rq, x, y.
    local ds_dx is EH_SurfaceDx(rp, rq, x, y).
    local dm_dx is self:k_fixed * (1 - x*y) * (1 - 3*x*y).

    return ds_dx - dm_dx.
}

// y gradient of alternative inequality
function EH_FixedTrackDy {
    parameter self, rp, rq, x, y.
    local ds_dy is EH_SurfaceDy(rp, rq, x, y).
    local dm_dy is -2 * self:k_fixed * x^2 * (1 - x*y).

    return ds_dy - dm_dy.
}

// Checks if the original inequality with variable track width has valid y
// solutions with rp & rq
// @param rp the rp value to test - co-prime to rq
// @param rq the rq value to test - co-prime to rp
// @return the found eccentricity limits, if any
function EH_FreeTrackLimits {
    parameter self, rp, rq.
    // find the lower limit on eccentricity by increasing from bottom
    local e_min is FindLimit(
        {parameter x, y. return EH_FreeTrack(self, rp, rq, x, y).},
        {parameter x, y. return EH_FreeTrackDx(self, rp, rq, x, y).},
        {parameter x, y. return EH_FreeTrackDy(self, rp, rq, x, y).},
        0
    ).

    // none only returned if continuous < 0 from top to bottom: no solutions
    if e_min = NONE {
        return NONE.
    }

    // find the upper limit on eccentricity by descending from top
    local e_max is FindLimit(
        {parameter x, y. return EH_FreeTrack(self, rp, rq, x, y).},
        {parameter x, y. return EH_FreeTrackDx(self, rp, rq, x, y).},
        {parameter x, y. return EH_FreeTrackDy(self, rp, rq, x, y).},
        1
    ).

    // max < min means no continuous band
    if e_max = NONE or e_max < e_min{
        return NONE.
    }

    return list(e_min, e_max).
}

// Solve variant equation with fixed track width. Will only be smaller
// than original above best altitude so validates solution in that
// domain.
//
// @param solution a soution found for the free track problem.
// @return the found eccentricity limits, if any
function EH_ApplyFixedTrackLimits {
    parameter self, solution.

    local rp is solution:rp.
    local rq is solution:rq.
    if solution:sma*(solution:e_min + 1) > self:alt_best + self:body_radius {
        // find the lower limit on eccentricity by increasing from bottom
        local e_min is FindLimit(
            {parameter x, y. return EH_FixedTrack(self, rp, rq, x, y).},
            {parameter x, y. return EH_FixedTrackDx(self, rp, rq, x, y).},
            {parameter x, y. return EH_FixedTrackDy(self, rp, rq, x, y).},
            0
        ).

        // none only returned if continuous < 0 from top to bottom: no solutions
        if e_min = NONE or e_min > solution:e_max {
            return NONE.
        }
        set solution:e_min to max(solution:e_min, e_min).
    }

    if solution:sma*(solution:e_max + 1) > self:alt_best + self:body_radius {
        // find the upper limit on eccentricity by descending from top
        local e_max is FindLimit(
            {parameter x, y. return EH_FixedTrack(self, rp, rq, x, y).},
            {parameter x, y. return EH_FixedTrackDx(self, rp, rq, x, y).},
            {parameter x, y. return EH_FixedTrackDy(self, rp, rq, x, y).},
            1
        ).

        // max < min means no continuous band
        if e_max = NONE or e_max < solution:e_min {
            return NONE.
        }
        set solution:e_max to min(solution:e_max, e_max).
    }

    return solution.
}

// Checks for hard limits on eccentricity with current sma.
// These are:
//     - required altitude over poles for scanning,
//     - minimum altitude at periapsis to stay safe
//     - maximum altitude at apoapsis to stay within scanning range and SOI
// @return the calculated hard limits on eccentricity, or None if they
// cannot be met
function EH_ApplyHardLimits {
    parameter self, solution.

    local a is solution:sma.
    local rr is self:body_radius.
    
    // altitude over poles above min altitude
    local e_pole is sqrt(1 - (rr + self:scanner:alt_min) / a).
    // apoapsis within max altitude
    local e_apo is (rr + self:scanner:alt_max) / a - 1.

    set solution:e_max to min(
        solution:e_max,
        min(e_pole, e_apo)
    ).
    if solution:e_min > solution:e_max {
        return NONE.
    }
    return solution.
}
// ################################ FUNCTIONS #################################

// get the greatest common denominator of a and b. Used to check for coprimes.
function gcd {
    parameter a, b.
    if a = 0 {
        return b.
    }
    if b = 0 {
        return a.
    }

    return gcd(b, mod(a, b)).
}

// Scale the fov to the target body as done inside ScanSat.
function scaleFov {
    parameter target_body, scanner.
    local fov is scanner:fov.
    local alt_best is scanner:alt_best.

    if target_body:radius < kerbin:radius {
        set fov to fov * sqrt(kerbin:radius / target_body:radius).
    }

    if fov > FOV_MAX {
        set alt_best to alt_best * FOV_MAX / fov.
        set fov to FOV_MAX.
    }

    return list(fov, alt_best).
}

// Finds a root between the specified starting values.
//
// Uses a bisection algorithm to quickly find a root between the two values.
// If multiple roots exist any one of them may be found.
//
// Assumes that f(x0) < 0 and f(x1) > 0 so should find roots that are
// increasing when moving from x0 to x1.
// x0 does not need to be smaller than x1 for this reason. If there are no
// roots, it will more towards closest value to 0, but not guaranteed.
//
// @param fx the function to find roots in
// @param x0 value where f(x) < 0
// @param x1 value were f(x) > 0
// @return x coordinate of found root
function findRootBetween {
    parameter fx, x0, x1.

    local y0 is fx(x0).
    local y1 is fx(x1).

    until abs(x0 - x1) < TOLERANCE {
        local x is (x0 + x1) / 2.
        local y is fx(x).
        if y < 0 {
            set x0 to x.
            set y0 to y.
        } else {
            set x1 to x.
            set y1 to y.
        }
    }

    if y0*y1 >= 0 {
        if abs(y0) <= abs(y1) {
            return x0.
        }
        return x1.
    }
    return (x0 + x1) / 2.
}

// Finds a root fx(x) = 0 near x0. Only searches in the direction specified.
//
// Uses a modified version of the Newtonian root finding algorithm to limit
// step size and force searching in a set direction from the starting value.
//
// Limit to step size was required to prevent overshoot taking it too far past
// the intended root and out of the domain.
//
// @param fx the function to find a root in
// @param df_dx derivative of f(x)
// @param x0 starting value
// @param d direction to search in (-1 for down, 1 for up)
// @return the 'x' value of the root
function findRootNear {
    parameter fx, df_dx, x0, d.

    local y0 is fx(x0).
    local _x is -1.
    local x is x0.
    until abs(_x - x) < TOLERANCE {
        
        local y is fx(x).
        local dy is df_dx(x).

        if y*y0 < 0 { // signs differ or one is root. bisect to root is faster
            if y0 < y {
                return findRootBetween(fx, _x, x).
            }
            return findRootBetween(fx, x, _x).
        }

        local rr is choose y/dy if dy <> 0 else 2*TOLERANCE.
        
        if abs(rr) > MAX_STEP {
            set rr to d * MAX_STEP.
        } else if rr*d < 0 {
            set rr to -rr.
        }

        set _x to x.

        set x to x + rr.

        if x < 0 or x > 1 {
            return NONE.
        } 
    }

    return x.
}

// Find the root with the largest y value that can be reached from the passed
// side of the domain (BOTTOM -> y = 0, TOP -> y = 1) and return its y value.
// If there are continuous roots between, None is returned instead.
//
// @param fxy The function to find roots in
// @param df_dx dirivative in x direction
// @param df_dy derivative in y direction
// @param side side of plot to start on (BOTTOM 0 or TOP 1)
// @return the y value of the root reached, or None if continuum
function findLimit {
    parameter fxy, df_dx, df_dy, side.

    local d is 1 - 2*side.

    // set starting position to corners of domain.
    local x is 1 - side.
    local y is side.


    if fxy(x, y) > 0 {  // starting above zero, move x to root
        set x to findRootNear(
            {parameter _x. return fxy(_x, y).},
            {parameter _x. return df_dx(_x, y).},
            x,
            -d
        ).


        if x = NONE {  // no root in x found, return current y.
            return y.
        }
    } else {  // starting below zero, move y to root
        set y to findRootBetween(
            {parameter _y. return fxy(x, _y).},
            y,
            1 - y
        ).

        
        // root y increases outside domain, y is limit.
        local dx is d * df_dx(x, y).

        if dx < 0 {
            return y.
        }
    }
    local x0 is x.
    local x1 is side.

    // repeat until x known to tolerance.
    until abs(x0 - x1) < TOLERANCE {
        local dx is d * df_dx(x, y).
        // decide which x needs updating to match root.
        if dx > 0 {
            set x0 to x.
            set x1 to findRootBetween(
                {parameter _x. return fxy(_x, y).},
                x0,
                x1
            ).
        } else {
            set x1 to x.
            set x0 to findRootBetween(
                {parameter _x. return fxy(_x, y).},
                x1,
                x0
            ).
        }

        // bisect x roots
        set x to (x0 + x1) / 2.
        // find value at new position (usually negative)
        local z is fxy(x, y).
        // move y to root at new x
        set y to findRootNear(
            {parameter _y. return fxy(x, _y).},
            {parameter _y. return df_dy(x, _y).},
            y,
            choose d if z < 0 else -d
        ).
        

        // y left domain
        if y = NONE {
            // if z was positive, left the way it entered so return that value
            return choose NONE if z < 0 else side.
        }
    }
    return y.
}


function freeTrackLimits {
    parameter eq_holders, rp, rq.

    local e_min is 0.
    local e_max is 1.
    for eq in eq_holders {
        if eq:alt_best > eq:sma_min - eq:body_radius {
            local e_lims is EH_FreeTrackLimits(eq, rp, rq).
            if e_lims = NONE {
                return NONE.
            }
            set e_min to max(e_min, e_lims[0]).
            set e_max to min(e_max, e_lims[1]).
            if e_min > e_max {
                return NONE.
            }
        }
    }
    return list(e_min, e_max).
}

function applyFixedTrackLimits {
    parameter eq_holders, solution.
    for eq in eq_holders {
        set solution to EH_ApplyFixedTrackLimits(eq, solution).
        if solution = NONE {
            return NONE.
        }
    }
    return solution.
}

function applyHardLimits {
    parameter eq_holders, solution.
    // altitude at periapsis
    local rr is eq_holders[0]:body_radius.
    local e_safe is 1 - (rr + eq_holders[0]:alt_safe) / solution:sma.
    // apoapsis within SOI
    local e_soi is eq_holders[0]:body:soiRadius / solution:sma - 1.

    set solution:e_max to min(solution:e_max, min(e_soi, e_safe)).

    for eq in eq_holders {
        set solution to EH_ApplyHardLimits(eq, solution).
        if solution = NONE {
            return NONE.
        }
    }
    return solution.
}

function findFastestOrbits {
    parameter target_body, alt_safe, ship_scanners.
    parameter update_col is -1.
    parameter update_row is -1.

    set DEBUG_ROW to choose update_row+1 if update_row <> -1 else -1.

    local eq_holders is list().
    local sma_max is target_body:soiRadius.
    local sma_min is target_body:radius.
    for scanner in ship_scanners {
        local eq is NewEquationHolder(target_body, alt_safe, scanner).
        eq_holders:add(eq).
        set sma_max to min(sma_max, eq:sma_max).
        set sma_min to max(sma_min, eq:sma_min).
    }

    local r_max is (eq_holders[0]:sma_sync / sma_max) ^ (3/2).
    local r_min is (eq_holders[0]:sma_sync / sma_min) ^ (3/2).

    local valid is list().
    local rp is 0.
    until valid:length > 0 {
        set rp to rp+1.
        local solutions is list().
        local rq is ceiling(rp*r_max) - 1.
        local end is floor(rp*r_min).
        until rq >= end {
            set rq to rq+1.
            if update_col <> -1 {
                print ("" + rp + "/" + rq):padright(terminal:width-update_col) at (update_col, update_row).
            }
            if gcd(rp, rq) = 1 {
                local e_limits is freeTrackLimits(eq_holders, rp, rq).
                if e_limits <> NONE {
                    local a is eq_holders[0]:getSma(rp, rq).
                    solutions:add(
                        NewSolutionParams(
                            rp,
                            rq,
                            a,
                            e_limits[0],
                            e_limits[1]
                        )
                    ).
                } else if solutions:length > 0 {
                    break.
                }
            }
        }

        set DEBUG_FLAG to true.
        local itr is solutions:reverseIterator.
        until not itr:next {
            local solution is itr:value.
            if update_col <> -1 {
                print ("validating " + solution:rp + "/" + solution:rq):padright(terminal:width-update_col) at (update_col, update_row).
            }
            set DEBUG_FLAG to true.
    
            set solution to applyHardLimits(eq_holders, solution).
            if solution <> NONE {
        
                set solution to applyFixedTrackLimits(eq_holders, solution).
                if solution = NONE {
                    break.
                }
        
                valid:add(solution).
            }
        }
        set DEBUG_FLAG to false.
    }

    return valid.
}


local DEBUG_FLAG is false.
function debug {
    parameter message.
    if not DEBUG_FLAG {
        return.
    } 


    print message:padRight(terminal:width) at (0, DEBUG_ROW).
    terminal:input:getChar().
    print "":padRight(terminal:width) at (0, DEBUG_ROW).
}

// ################################# OBJECTS ##################################
// create instances for all scanners included in scansat
global SCANNERS is lex(
    "ms-1", NewScanner(3, 20_000, 70_000, 250_000), 
    "ms-2a", NewScanner(4, 100_000, 500_000, 750_000), 
    "ms-r", NewScanner(1.5, 70_000, 300_000, 400_000), 
    "r-3b", NewScanner(1.5, 5_000, 70_000, 250_000), 
    "r-eo-1", NewScanner(3.5, 50_000, 100_000, 500_000), 
    "sar-c", NewScanner(3, 500_000, 700_000, 750_000), 
    "sar-l", NewScanner(4, 250_000, 500_000, 1_000_000), 
    "sar-x", NewScanner(1.5, 70_000, 250_000, 500_000), 
    "scan-r", NewScanner(1, 20_000, 70_000, 250_000), 
    "scan-r2", NewScanner(2.5, 70_000, 250_000, 500_000), 
    "scan-rx", NewScanner(3, 100_000, 500_000, 750_000), 
    "vs-1", NewScanner(1.5, 20_000, 70_000, 250_000), 
    "vs-11", NewScanner(4, 100_000, 200_000, 1_000_000), 
    "vs-3", NewScanner(2.5, 70_000, 350_000, 500_000) 
).