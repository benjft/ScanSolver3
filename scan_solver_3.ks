local BEEP is char(7).

function input {
    parameter col, row, value is "", decal is ">>".
    local ti is terminal:input.
    until false {
        local dispStr is decal + value.
        print dispStr:padRight(terminal:width) at (col,row).
        local c is ti:getChar.
        if c = ti:backspace {
            set value to value:subString(0, max(0, value:length-1)).
        } else if c = ti:enter {
            return value.
        } else if unChar(c) > 31 and unChar(c) < 127 {
            set value to value + c.
        }
    }
}

function scroll_view {
    parameter col, row, strings.
    local ti is terminal:input.

    local i is 0.
    local n is strings:length.

    until false {
        local dispStr is strings[i].
        local idxStr is "" + (i+1) + " of " + n + ": ".

        print (idxStr + dispStr):padRight(terminal:width) at (col, row).
        local c is ti:getChar.

        if c = ti:upCursorOne {
            set i to mod(i-1, n).
        } else if c = ti:downCursorOne {
            set i to mod(i+1, n).
        } else if c = ti:rightCursorOne {
            set i to mod(i - 10, n).
        } else if c = ti:leftCursorOne {
            set i to mod(i + 10, n).
        } else if c = ti:enter {
            break.
        }
        until i >= 0 {
            set i to i + n.
        }
    }
}

function firstMissing {
    parameter list0, list1.

    for item in list0 {
        if list1:find(item) = -1 {
            return item.
        }
    }
    return 0.
}

runOncePath("./lib_scan_solver_3.ks").
function main {
    set terminal:width to 80.
    set terminal:height to 40.

    clearScreen.

    local row is 0.
    print "Input Target Body:" at (0, row).
    local target_body is input(0, row+1).
    until bodyExists(target_body) {
        print ("cannot find " + target_body + ". Check spelling!"):padRight(terminal:width) at (0, row+1).
        wait until terminal:input:getchar() = terminal:input:enter.
        set target_body to input(0, row+1).
    }
    set target_body to body(target_body). 

    set row to row + 3.
    print "Input Minimum Safe Altitude: (meters)" at (0, row).
    local alt_safe is input(0, row+1).
    until alt_safe:toNumber(-1) >= 0 {
        print "Positive number expected":padRight(terminal:width) at (0, row+1).
        wait until terminal:input:getchar() = terminal:input:enter.
        set alt_safe to input(0, row+1).
    }
    set alt_safe to alt_safe:toNumber.

    set row to row + 3.
    print "Input Scanner Names (short, space separated)" at (0, row).
    local name is input(0, row+1).
    local names is name:split(" ").
    local nf is firstMissing(names, SCANNERS:keys).
    until nf = 0 {
        print "Name '" + nf + "' not found!.":padRight(terminal:width) at (0, row+1).
        wait until terminal:input:getchar() = terminal:input:enter.
        set name to input(0, row+1, name).
        set names to name:split(" ").
        set nf to firstMissing(names, SCANNERS:keys).
    }

    local scannerList is list().
    for name in names {
        scannerList:add(SCANNERS[name]).
    }

    local old_ipu is config:ipu.
    set config:ipu to 2000.
    set row to row + 3.
    print "Working...":padright(terminal:width) at (0, row).
    local t is time.
    local solutions is findFastestOrbits(target_body, alt_safe, scannerList, "Working... ":length, row).
    set t to time - t.
    set config:ipu to old_ipu.

    print BEEP.
    print ("completed calculation in " + round(t:seconds, 1) + " seconds"):padright(terminal:width) at (0, row+2).
    
    local strings is list().
    for solution in solutions {
        local rp is "" + solution:rp.
        local rq is "" + solution:rq.
        local e_min is "" + CEILING(solution:e_min, 5).
        local e_max is "" + floor(solution:e_max, 5).
        local a is "" + round(solution:sma, 2).

        local str is "(" + rp + "/" + rq + ")" + " a=" + a + "m  e=" + e_min +
            " to " + e_max.
        strings:add(str).
    }

    print "Use arrow keys to scroll. Press enter to end.":padright(terminal:width) at (0, row+1).
    scroll_view(0, row, strings).

    return row+4.
}

LOCAL oldIPU IS CONFIG:IPU.
SET CONFIG:IPU TO 2000. //2000 is the max allowed IPU by the kOS source code
local _row is main().
print "start again? (y/n)" at (0, _row+1).
until input(0, _row+2):startsWith("n") {
    set _row to main().
    print "start again? (y/n)" at (0, _row+1).
}
clearScreen.
SET CONFIG:IPU TO oldIPU.