#!/usr/bin/env python
r"""
Enhanced solving of easier exercises from mathematical analysis
(aimed at the difficulty level of freshmen students at the Faculty of
Mathematics and Physics of the Charles University in Prague, Czech Republic)
AUTHORS:

- Karel Ha (2012): initial version (with kind help of supervisor Robert Samal)
- Robert Samal (2012): suggestion for message "For usage tips run
  help_sage_alpha()" 


"""
#*****************************************************************************
#      Copyright (C) 2012 Karel Ha <mathemage@gmail.com>
#
# Distributed under the terms of the GNU General Public License (GPL)
# as published by the Free Software Foundation; either version 2 of
# the License, or (at your option) any later version.
#                 http://www.gnu.org/licenses/
#*****************************************************************************

# WELCOME MESSAGE
print r"""SageAlpha, version 1.0.0
http://ha.matfyz.cz/sagealpha

********************************************************************
      Copyright (C) 2012 Karel Ha <mathemage@gmail.com>

 Distributed under the terms of the GNU General Public License (GPL)
 as published by the Free Software Foundation; either version 2 of
 the License, or (at your option) any later version.
                 http://www.gnu.org/licenses/
********************************************************************

For usage tips run help_sage_alpha()
"""

###################################
# SAMPLE FUNCTIONS
myfunc0 = 4.2
myfunc1 = 1/4*(x - 4)*(pi + 4*x)*x
myfunc2(y) = (2*(y+.5)^2 + 3*(y+.5) + 1)/(y+.5)
myfunc3 = arctan(3*x) + 6 - pi/2
myfunc4 = 1/(x-2) + 3
myfunc5(y) = y^2 + 3*y
myfunc6 = x+sin(x)
myfunc7 = sin(x/pi)
myfunc8 = log(x) + sin(x)
myfunc9 = sin(x^2)
myfunc10 = x*(x-2)/((x-2)*(x+10)*(x-3)*log(abs(x)))         # too brutal :-)
myfunc11(x) = x*abs(log(abs(x)))
###################################

def help_sage_alpha():
    r"""
    This command shows a brief list of function available in SageAlpha module
    """
    print r"""  For given function "f" following commands are available:
asymptotes()                - "tangent" to a curve at infinity/limitly at a certain point
curvature()                 - intervals where "f" is convex / concave
inflection_points()         - points with zero 2nd derivative
investigations()            - dictionary with various findings about "f" (see other commands in this section)
monotonicity()              - intervals where "f" is rising / declining
print_investigations()      - output to the console / notebook, in text form
stationary_points()         - points with zero 1st derivative
*tex_investigations()       - output to the pdf document, in LaTeX-ed form
texie()                     - same as "tex_investigations()", with "okular" as pdf viewer & LaTeX's source code printed out
undefined()                 - x-coordinates for undefined points
________________________________________________________________
  Extra functionality is provided with these auxiliary commands:
all_print_investigations()  - consecutive test-runs of sample functions
*abs_args()                 - get recursively operands of "abs" subfunctions
*abs_eval()                 - recursively apply the given action to formulae acquired by elimination of "abs" occurrence
extract_args()              - get recursively operands of the given unary operator (in subformulae)
extrapolate_period()        - integer period of the given formula
tangent()                   - render tangent to "f" at the given point in the plot
tex_list()                  - LaTeX's code for the set of given formulae

* under construction

For further details on a command, run "command?" """

#######################################################
# test all sample functions with "print_investigations"
def all_print_investigations(upto=11, prefix="myfunc"):
    myfuncs = [prefix+str(i) for i in range(upto+1)]
    for i in range(upto+1): print_investigations(eval(myfuncs[i]))
#######################################################

# when limits are not real numbers
bad_limits = [SR('und'), SR('ind'), -oo, oo]

def undefined(self):
    r"""
    This function returns list of x-coordinates for undefined points of ``self``
    (i.e. where function is out of its domain).

    INPUT:
        
    - ``self`` - the tested real-valued function ``f``
        This should be a function from ``R`` numbers to ``R`` numbers

    OUTPUT:

    a list  -- x-coordinates of undefined points

    EXAMPLES:

    This example illustrates default usage::

        sage: undefined(1/(x - 2) + 3)
        [2]

    NOTES:

    This function first finds points lying out of function's domain. Initial
    version works for simpler (e. g. rational, log()) functions through finding
    null-points of denominators and of expressions in log() arguments.

    TODO:

    - further causes of ``undefinition``: negative arguments of square roots
      (and their fourth, sixth root analogies...), arcsin & arccos arguments...  

    AUTHORS:

    - Karel Ha (2012-04-07): initial version
    """
    if self.is_real():              # a constant function
        return []

    # "f(x) = sth(x)" form instead of "f = sth(x)"
    x = var('x')
    assume(x, 'real')
    fnct(x) = self.subs({self.default_variable():x})

    eqns = solve(fnct.simplify_full().denominator() == 0, x,
            to_poly_solve=True)
    tmp = [eqn.right() for eqn in eqns if eqn.left() == x]

    eqns2 = solve(fnct.simplify_full().denominator() == 0, x,
            to_poly_solve="force")
    tmp += [eqn.right() for eqn in eqns2 if eqn.left() == x]

    undefs = list(set(tmp))
    undefs.sort()
    forget(x, 'real')
    return undefs

def stationary_points(self):
    r"""
    This function returns list of tuples with coordinates of stationary points
    belonging to a real-valued function (of 1 variable).

    INPUT:
        
    - ``self`` - the tested real-valued function ``f``
        This should be a function from ``R`` numbers to ``R`` numbers

    OUTPUT:

    - a list of tuples -- pairs "(x, f(x))" for ``x`` as one of
        stationary points

    EXAMPLES:

    This example illustrates default usage::

        sage: g = (y-1)*(y+3)
        sage: list = stationary_points(g); list
        [(-1, -4)]

    NOTES:

    This function first finds domain points where 1st derivative equals zero.
    Consequently, it creates a list of respective pairs of aforementioned
    points and their functional values. This is done while checking the domain
    point is of type "(x : real_value)".    

    TODO:

    - stationary points for x-values where 1st derivative is not defined

    AUTHORS:

    - Karel Ha (2012-03-30): initial version
    """
    # "f(x) = sth(x)" form instead of "f = sth(x)"
    # "fnct" will now use default variable "x"
    if self.is_real():              # CONSTANT FUNCTION - TEMPORARY SOLUTION!!
        return [(-oo,oo)]
    else:
        #x = var('x')
        x = self.default_variable()
        assume(x, 'real')
        fnct(x) = self.subs({self.default_variable():x})

    try:
        tmp = solve([fnct.diff(1) == 0], x, to_poly_solve=True)
    except:
        tmp = solve([fnct.diff(1) == 0], x)

    tmp += solve([fnct.diff(1) == 0], x, to_poly_solve="force")
    x_coords = list(set(tmp))

    # non-periodic points
    stats = [(pt.right(), fnct().subs({x:pt.right()}).simplify_full())\
            for pt in x_coords\

            # no general meaningless solutions
            # +0*pi ... hack for converting to symbolic expressions (numbers
            # don't have args() method)
            if pt.left() == x and (pt.right()+0*pi).args() == () and\

            # real-valued points
            pt.right().n() != NaN and pt.right().n().is_real()]

    # periodic solutions: simplified functional values
    if var('k') not in self.args():         # dummy integer variable
        k = var('k')
    else:
        k = var('n')
    assume(k, 'integer')

    for pt in x_coords:
        if pt.left() == x and len(pt.right().args()) > 0 and x not in \
        pt.right().args() and pt.right()(0).imag().simplify_full() == 0:
            ptright = pt.right().subs({ pt.right().args()[0] : k })
            stats += [(ptright,
                fnct().subs({x:ptright}).simplify_full())]

    forget(k, 'integer')
    stats.sort()
    forget(x, 'real')
    return stats

def inflection_points(self):
    r"""
    This function returns list of tuples with coordinates of inflection points
    belonging to a real-valued function (of 1 variable).

    INPUT:
        
    - ``self`` - the tested real-valued function ``f``
        This should be a function from ``R`` numbers to ``R`` numbers

    OUTPUT:

    a list of tuples -- pairs "(x, f(x))" for ``x`` as one of
        inflection points

    EXAMPLES:

    This example illustrates default usage::

        sage: g = (y-1)*(y+3)*y
        sage: ls = inflection_points(g); ls    
        [(-2/3, 70/27)]

    ::

        sage: inflection_points(cos(x))
        [(1/2*pi + 2*pi*z67, 0), (-1/2*pi + 2*pi*z65, 0)]

    NOTES:

    This function first finds domain points where 2nd derivative equals zero.
    Consequently, it creates a list of respective pairs of aforementioned
    points and their functional values. This is done while checking the domain
    point is of type "(x : real_value)".    
    
    TODO:

    - inflection points for x-values where 2nd derivative is not defined

    AUTHORS:

    - Karel Ha (2012-03-31): initial version
    """
    # "f(x) = sth(x)" form instead of "f = sth(x)"
    x = var('x')
    assume(x, 'real')
    # "fnct" will now use default variable "x"
    if self.is_real():              # CONSTANT FUNCTION - TEMPORARY SOLUTION!!
        return [(-oo,oo)]
    else:
        fnct(x) = self.subs({self.default_variable():x})

    try:
        tmp = solve([fnct.diff(2) == 0], x, to_poly_solve=True)
    except:
        tmp = solve([fnct.diff(2) == 0], x)

    tmp += solve([fnct.diff(2) == 0], x, to_poly_solve="force")
    x_coords = list(set(tmp))

    # non-periodic points
    infls = [(pt.right(), fnct().subs({x:pt.right()}).simplify_full())\
            for pt in x_coords\

            # no general meaningless solutions
            # +0*pi ... hack for converting to symbolic expressions (numbers
            # don't have args() method)
            if pt.left() == x and (pt.right()+0*pi).args() == () and\

            # real-valued points
            pt.right().n() != NaN and pt.right().n().is_real()]

    # periodic solutions: simplified functional values
    if var('k') not in self.args():         # dummy integer variable
        k = var('k')
    else:
        k = var('n')
    assume(k, 'integer')

    for pt in x_coords:
        if pt.left() == x and len(pt.right().args()) > 0 and x not in \
        pt.right().args() and pt.right()(0).imag().simplify_full() == 0:
            ptright = pt.right().subs({ pt.right().args()[0] : k })
            infls += [(ptright,
                fnct().subs({x:ptright}).simplify_full())]

    forget(k, 'integer')
    infls.sort()
    forget(x, 'real')
    return infls

def asymptotes(self, astype="all"):
    r"""
    This function returns the function's list of asymptotes (horizontal,
    vertical and oblique).

    INPUT:
        
    - ``self`` - the input real-valued function
        This should be a function from ``R`` numbers to ``R`` numbers, even
        constant function should at least be in form x `|-->` 42
    - ``astype`` - asymptote's type (default: "all") - "all",
      "oblique", "vertical" or "horizontal"

    OUTPUT:

    a list -- functions of asymptotes

    EXAMPLES:

    This example illustrates default (i.e. showing all) usage::

        sage: asymptotes(myfunc2)
        [x |--> 2.0*x + 4.0, y |--> -1/2]

    Showing only vertical asymptotes::

        sage: asymptotes(myfunc2, astype="vertical")
        [y |--> -1/2]

    Showing only oblique asymptotes::

        sage: asymptotes(myfunc2, astype="oblique")
        [x |--> 2.0*x + 4.0]

    Showing only horizontal asymptotes::

        sage: asymptotes(myfunc3, astype="horizontal")
        [x |--> -pi + 6, x |--> 6]

    TODO:

    - verify proper functionality of periodical (vertical) asymptotes - e. g.
      in tan(x)

    AUTHORS:

    - Karel Ha (2012-03-31): initial version
    """
    x = var('x')
    if astype in ["horizontal", "all"] and self() in RR:     # constant function
        f(x) = self()
        return [f]

    # "f(x) = sth(x)" form instead of "f = sth(x)"
    assume(x, 'real')
    fnct(x) = self.subs({self.default_variable():x})

    results = []

    if astype in ["horizontal", "all"]:
        # horizontal asymptote -> -oo in form of y(x) = c
        c = lim(fnct, x = -oo)
        if c not in bad_limits:
            results.append(c + 0*x)

        # horizontal asymptote -> oo in form of y(x) = c
        c = lim(fnct, x = oo)
        if c not in bad_limits:
            results.append(c + 0*x)

    if astype in ["oblique", "all"]:
        # oblique asymptote -> -oo in form of y(x) = ax+b 
        a = lim(fnct/x, x = -oo)
        # proper lim & non-constant asymptote
        if a not in bad_limits and a != 0:
            b = lim(fnct - a*x, x = -oo)
            if b not in bad_limits:
                asy(x) = a*x + b
                results.append(asy)

        # oblique asymptote -> oo in form of y(x) = ax+b 
        a = lim(fnct/x, x = oo)
        # proper lim & non-constant asymptote
        if a not in bad_limits and a != 0:
            b = lim(fnct - a*x, x = oo)
            if b not in bad_limits:
                asy(x) = a*x + b
                results.append(asy)

    # get rid of duplicities
    results = list(set(results))

    if astype in ["vertical", "all"]:
        # vertical asymptotes: undefined points & one-sided limits leads to +- infty
        undefs = undefined(self)
        for x_coord in undefs:
            if lim(fnct, x = x_coord, dir='-') in [-oo,oo]\
               or lim(fnct, x = x_coord, dir='+') in [-oo,oo]:
                y = var('y')
                f(y) = x_coord
                results += [f]

    forget(x, 'real')
    return results

def monotonicity(self, type="inc", hints=[]):
    r"""
    This function returns list of conditions on ``x`` for non-decreasingness or
    non-increasingness of ``self``.

    INPUT:
        
    - ``self`` - the input real-valued function ``f``
        This should be a function from ``R`` numbers to ``R``
    - ``type`` - string "inc" or "dec" (default: "inc") 
        If "inc", search for strictly increasing parts, else for strictly
        decreasing parts.
    - ``hints`` - list (default: "[]") 
        Suggestions from the user for end-points of tested intervals (including
        infinities).

    OUTPUT:

    a list - intervals (tuples of endpoints) ``x`` for given type of
    monotonicity

    EXAMPLES:

    This example illustrates default (i.e. increasing monotonicity) usage::

        sage: monotonicity(x^4-3*x^2+2)      
        [(-1/2*sqrt(2)*sqrt(3), 0), (1/2*sqrt(2)*sqrt(3), +Infinity)]

    This example illustrates usage of decreasing monotonicity::

        sage: monotonicity(x+1/x, type="dec") 
        [(-1, 0), (0, 1)]

    This example displays behavior for incorrect input::

        sage: monotonicity(x+1/x, type="anti")   
        (...)
        TypeError: "anti" is not a valid type of monotonicity.

    NOTES:

    Firstly, points splitting real numbers into intervals are acquired (either
    from ``hints`` or manually through ``stationary_points()``). From each
    interval one ``test_point`` is chosen at which derivative is evaluated,
    thus, its increasingness/decreasingness is recognized.

    TODO:

    - either improve splitting points by adding ``stationary_points`` where 1st
      derivative is undefined
    - or find better way to find these intervals (preferably using "solve" that
      will find intervals of positive/negative 1st derivative)

    AUTHORS:

    - Karel Ha (2012-04-02): initial version
    """
    if type not in ["inc", "dec"]:
        raise TypeError ('"'+type+"\" is not a valid type of monotonicity.")

    # endpoints of intervals
    if hints == []:
        end_pts = [-oo]
        end_pts.extend([pt[0] for pt in stationary_points(self) \
            if pt[0].args()==() ])
        end_pts.extend(undefined(self))
        end_pts += [oo]
    else:
        end_pts = hints

    no_infties = set(end_pts)
    no_infties.discard(-oo)
    no_infties.discard(oo)

    tmp_ls = sorted(list(no_infties))
    if -oo in end_pts: tmp_ls = [-oo] + tmp_ls
    if oo in end_pts: tmp_ls = tmp_ls + [oo]
    end_pts = tmp_ls

    if self().args() == ():                     # constant function
        # ...is neither strictly increasing nor decreasing
        return [(-oo,+oo)]
    else:                                       # non-constant function
        # tested intervals of monotonicity
        ints = [ (end_pts[i-1],end_pts[i]) for i in range(1,len(end_pts))]

        # representants for tested intervals of monotonicity
        test_pts = []
        for it in ints:
            if [it[0], it[1]] == [-oo, oo]:
                test_pts.append(0)
            elif it[0] == -oo:
                test_pts.append(it[1]-1)
            elif it[1] == oo:
                test_pts.append(it[0]+1)
            else:
                test_pts.append((it[0]+it[1])/2)

        # form of derivative's function with 'x' variable
        x = var('x')
        Dself(x) = self.diff(1)(x)

        if type == "inc":
            op = ">"
        elif type == "dec":
            op = "<"
        condition = "Dself().subs({x:test_pts[i]})"+op+"0"
        return [ints[i] for i in range(len(ints)) if eval(condition)]

def curvature(self, type="convex", hints=[]):
    r"""
    This function returns list of conditions on ``x`` for convexity or concavity
    of ``self``.

    INPUT:
        
    - ``self`` - the input real-valued function ``f``
        This should be a function from ``R`` numbers to ``R``
    - ``type`` - string "convex" or "concave" (default: "convex") 
        If "convex", search for convex parts, else for concave parts.
    - ``hints`` - list (default: "[]") 
        Suggestions from the user for end-points of tested intervals (including
        infinities).

    OUTPUT:

    a list - intervals (tuples of endpoints) ``x`` for given type of curvature

    EXAMPLES:

    This example illustrates default (i.e. convexity) usage::

        sage: curvature(x*(x-1)*(x+2))
        [(-1/3, +Infinity)]

    This example illustrates usage of concavity::

        sage: curvature(x*(x-1)*(x+2), type="concave") 
        [(-Infinity, -1/3)]

    This example displays behavior for incorrect input::

        sage: curvature(x*(x-1)*(x+2), type="cave")
        (...)
        TypeError: "cave" is not a valid type of curvature.

    NOTES:

    Firstly, points splitting real numbers into intervals are acquired (either
    from ``hints`` or manually through ``inflection_points()``). From each
    interval one ``test_point`` is chosen at which 2nd derivative is evaluated,
    thus, its convexity/concavity is recognized.

    TODO:

    - either improve splitting points by adding ``inflection_points`` where 2nd
      derivative is undefined
    - or find better way to find these intervals (preferably using "solve" that
      will find intervals of positive/negative 2nd derivative)

    AUTHORS:

    - Karel Ha (2012-04-02): initial version
    """
    if type not in ["convex", "concave"]:
        raise TypeError ('"'+type+"\" is not a valid type of curvature.")

    # endpoints of intervals
    if hints == []:
        end_pts = [-oo]
        end_pts.extend([pt[0] for pt in inflection_points(self)])
        end_pts.extend(undefined(self))
        end_pts += [oo]
    else:
        end_pts = hints

    no_infties = set(end_pts)
    no_infties.discard(-oo)
    no_infties.discard(oo)

    tmp_ls = sorted(list(no_infties))
    if -oo in end_pts: tmp_ls = [-oo] + tmp_ls
    if oo in end_pts: tmp_ls = tmp_ls + [oo]
    end_pts = tmp_ls

    if self().args() == ():                     # constant function
        # ...is neither strictly convex nor concave
        return [(-oo,+oo)]
    else:                                       # non-constant function
        # tested intervals of monotonicity
        ints = [ (end_pts[i-1],end_pts[i]) for i in range(1,len(end_pts))]

        # representants for tested intervals of curvature
        test_pts = []
        for it in ints:
            if [it[0], it[1]] == [-oo, oo]:
                test_pts.append(0)
            elif it[0] == -oo:
                test_pts.append(it[1]-1)
            elif it[1] == oo:
                test_pts.append(it[0]+1)
            else:
                test_pts.append((it[0]+it[1])/2)

        # form of 2nd derivative's function with 'x' variable
        x = var('x')
        Dself2(x) = self.diff(2)(x)

        if type == "convex":
            op = ">"
        elif type == "concave":
            op = "<"
        condition = "Dself2().subs({x:test_pts[i]})"+op+"0"
        return [ints[i] for i in range(len(ints)) if eval(condition)]

def tangent(self, point):
    r"""
    This function returns a tangent line of ``self`` at ``point``.

    INPUT:
        
    - ``self`` - the input real-valued function ``f``
        This should be a function from ``R`` numbers to ``R``
    - ``point`` - x-coordinate for point of touch

    OUTPUT:

    a function - linear function (a*x+b) of tangent if it exists; otherwise
    ``[]``

    EXAMPLES:

    This example illustrates default usage::

        sage: tangent(tan(x), 42*pi)
        x \|--> -42*pi + x

    AUTHORS:

    - Karel Ha (2012-04-13): initial version
    """
    fnct(x) = self.subs({self.default_variable():x})
    if point in undefined(fnct.diff()):
        return []
    else:
        slope = fnct.diff().subs({x:point})
        shift = fnct().subs({x:point})
        result(x) = slope*(x-point) + shift
        return result

def extrapolate_period(self):
    r"""
    This function struggles to find out the period of given formula.

    INPUT:
        
    - ``self`` - formula: expression with integer variable representing
      periodicity

    OUTPUT:

    a number/empty list - value expressing the period if successful; otherwise
    0

    EXAMPLES:

    This example illustrates default usage::

        sage: m = var('m')
        sage: assume(m, 'integer')
        
        sage: extrapolate_period(sqrt(-1/2*pi) + 2*pi*m)
        2*pi

    ::

        sage: extrapolate_period(sqrt(-1/2*pi + 2*pi*m))
        0

    AUTHORS:

    - Karel Ha (2012-04-15): initial version
    """

    x = self.default_variable()
    if self.is_polynomial(x) and self.degree(x) == 1:
        return self.subs({x:1}) - self.subs({x:0})
    else:
        return 0

def extract_args(self, operation, atomic=False):
    r"""
    This function extracts all arguments (in form of symbolic expressions) of
    unary ``operation()`` subfunctions of real-valued function ``self``

    INPUT:
        
    - ``self`` - symbolic expression of real-valued function
    - ``operand`` - string of unary real-valued operation
    - ``atomic`` - boolean (default: False)
        If true, return only expression without ``operation()`` subexpressions

    OUTPUT:

    a list - expressions of ``operation()`` subfunctions' arguments

    EXAMPLES:

    This example illustrates default usage::

        sage: extract_args(abs(x+abs(abs(abs(x/log(sqrt(log(log(x)))))))), abs)
        [x + abs(x/log(sqrt(log(log(x))))), x/log(sqrt(log(log(x))))]

    ::

        sage: extract_args(abs(x+abs(abs(abs(x/log(sqrt(log(log(x)))))))), log)
        [sqrt(log(log(x))), log(x), x]

    ::

        sage: extract_args(abs(x+abs(abs(abs(x/log(sqrt(log(log(x)))))))),
        arcsin)             
        []


    This example illustrates ``atomic`` switch::

        sage: extract_args(abs(x+abs(abs(abs(x/log(sqrt(log(log(x)))))))), abs,
        atomic=True)
        [x/log(sqrt(log(log(x))))]

    ::

        sage: extract_args(abs(x+abs(abs(abs(x/log(sqrt(log(log(x)))))))), log,
        atomic=True)
        [x]

    ::

        sage: extract_args(abs(x+abs(abs(abs(x/log(sqrt(log(log(x)))))))),
        sqrt, atomic=True)
        [log(log(x))]

    AUTHORS:

    - Karel Ha (2012-05-14): initial version
    """
    w0 = SR.wild()
    w1 = SR.wild()
    res = [f.operands()[0] for f in self.find(operation(w0))]

    subexprs = []
    for f in res:
        subexprs.extend( extract_args(f, operation, atomic) )

    if atomic:
        return [expr for expr in res+subexprs if not expr.has(operation(w1))]
    else:
        return res+subexprs

def abs_args(self, atomic=False):
    r"""
    This function extracts all arguments (in form of symbolic expressions) of
    ``abs()`` subfunctions of real-valued function ``self``

    INPUT:
        
    - ``self`` - symbolic expression of real-valued function
    - ``atomic`` - boolean (default: False)
        If true, return only expression without ``abs()`` subexpressions

    OUTPUT:

    a list - expressions of ``abs()`` subfunctions' arguments

    EXAMPLES:

    This example illustrates default usage::

        sage: c = abs(x+42)*x*abs(log(abs(x+abs(x-42+sin(x))+e^(abs(log(x)^3*cos(sin(x)))))))

    ::

        sage: abs_args(c)
        [x + 42, log(abs(x + e^(abs(log(x)^3*cos(sin(x)))) + abs(x + sin(x) - 42))),
        x + e^(abs(log(x)^3*cos(sin(x)))) + abs(x + sin(x) - 42), x + sin(x) - 42,
        log(x)^3*cos(sin(x))]

    This example illustrates ``atomic`` switch::

        sage: abs_args(c, atomic=True)
        [x + 42, x + sin(x) - 42, log(x)^3*cos(sin(x))]

    This example illustrates ability of merging multiple consecutive ``abs``
    expressions::

        sage: abs_args(abs(x+abs(abs(abs(x/log(x))))))
        [x + abs(x/log(x)), x/log(x)]

    AUTHORS:

    - Karel Ha (2012-04-22): initial version
    """
    return extract_args(self, abs, atomic)

def abs_eval(self, processor='print self.simplify_full()', assums=[]):
    r"""
    This function recursively eliminates ``abs()`` subfunctions. When done,
    ``processor`` is run for each interval/assumption where value of respective
    ``abs()`` is either non-negative or negative.

    INPUT:
        
    - ``self`` - symbolic expression of real-valued function
    - ``processor`` - Sage/Python function (default: print)
    - ``assums`` - list (default: empty list) assumptions from previous levels
      of recursion

    EXAMPLES:

    This example illustrates default usage::

        sage: abs_eval(abs(abs(x+2)-1), 'print self.simplify_full(); print
        self.simplify_full().diff()')
        [x + 2 >= 0, x + 1 >= 0]
        x |--> x + 1
        x |--> 1
        --------------------
        [x + 2 >= 0, x + 1 < 0]
        x |--> -x - 1
        x |--> -1
        --------------------
        [x + 2 < 0, x + 3 >= 0]
        x |--> x + 3
        x |--> 1
        --------------------
        [x + 2 < 0, x + 3 < 0]
        x |--> -x - 3
        x |--> -1
        --------------------

    TODO:

    - Not working - due to certain issues with inequalities in Sage 4.8. This
      problem is discussed & its solution can be viewed at:
      https://groups.google.com/forum/#!topic/sage-support/H-QcXkNCajM
    - Get rid of extra stuff in code (originally for the purpose of debugging)

    AUTHORS:

    - Karel Ha (2012-05-02): initial version - STILL UNDER CONSTRUCTION
    """
    #print assumptions()
    print "--------------------"
    
    abss = abs_args(self, atomic=True)
    if len(abss) == 0:
        print assumptions()
        #print self.simplify_full()
        exec(processor)
        print "--------------------"
    else:
        print "old assumptions()   ", assumptions()
        #new_assums = assums + [abss[0] >= 0]
        #assume(new_assums)
        print type(abss[0] > 0) 
        print abss[0] > 0
        assume( abss[0] > 0 )
        print "new assumptions()   ", assumptions()
        #print "new_assums          ", new_assums
        #print "assums              ", assums
        print "self                ", self
        print "self.simplify_full()", self.simplify_full()
        abs_eval(self.simplify_full(), processor, assums=[])
        #abs_eval(self.simplify_full(), processor, assums=new_assums)
        forget(abss[0] > 0)

        print "old assumptions() 2 ", assumptions()
        #new_assums = assums + [abss[0] < 0]
        #assume(new_assums)
        assume( abss[0] < 0 )
        print "new assumptions() 2 ", assumptions()
        print "self                ", self
        print "self.simplify_full()", self.simplify_full()
        abs_eval(self.simplify_full(), processor, assums=[])
        #abs_eval(self.simplify_full(), processor, assums=new_assums)
        forget(abss[0] < 0)

def investigations(self, left=-2*pi, right=2*pi, down=-15, up=15, tangents=True,
        autozoom=0):
    r"""
    This function returns a variety of information on given real-valued
    function, e. g. stationary/inflection points, extremes, asymptotes, its plot
    etc.

    INPUT:
        
    - ``self`` - the input real-valued function ``f``
        This should be a function from ``R`` numbers to ``R``
    - ``left`` - left margin of graph
    - ``right`` - right margin of graph
    - ``up`` - up margin of graph
    - ``tangents`` - bool (default: True); if True, display graph with tangent
      lines at the inflection points
    - ``autozoom`` - float (default: 1.1); if True, rescale the graph by the
      factor of ``autozoom`` in order to show as much important information
      (e. g. all reasonable stationary/inflection points, vertical asymptotes
      etc.)

    OUTPUT:

    a dictionary

    - keys: what was investigated
    - values: list of results from functions such as ``stationary_points``,
      ``inflection_points``, ``monotonicity``, ``period`` etc.

    TODO:

    - simplify complex formulae of points for rendering in the plot, see::

        investigations(x^5-x^3-x^2+x+1)

      due to troubles in:: 

        stationary_points(x^5-x^3-x^2+x+1)

    - Needs beta-testing for rendering of periodical (vertical) asymptotes
    - More comprehensible comments/documentation for ``autozoom`` feature
      (collecting positions of right, left, top, bottom borders; consequent
      rendering of whole plot)

    AUTHORS:

    - Karel Ha (2012-05-19): initial version
    - Robert Samal (2012-04-10): choice of color for tangent lines
    """
    outlist = []
    ###########################################################################
    # "F(x) = STH(x)" FORM INSTEAD OF "F = STH(x)" OR "F = STH(y)"
    x = var('x')
    if not self.is_real():                  # non-constant function
        _x_ = self.default_variable()
        fnct(x) = self.subs({_x_:x})
    else:                                   # constant function
        fnct(x) = (self+0*x)
    assume(x, 'real')
    fnct = fnct.full_simplify()
    ###########################################################################

    pic = Graphics()
    
    # coordinates for autozoom
    ideal_left = oo
    ideal_right = -oo
    ideal_down = oo
    ideal_up = -oo

    ###########################################################################
    # SET MARGINS (IDEAL_LEFT, IDEAL_UP...) FOR DISPLAYING GRAPH
    # ...AND GATHER CANDIDATES FOR THE PERIOD

    pers = []                                  # list of candidates for periods

    stats = stationary_points(fnct)
    outlist += [('stats', stats)]
    if stats not in [[], [(-oo,oo)]]:
        for pt in stats:
            if pt[0] in RR:                    # non-periodic solution
                pic += point(pt, rgbcolor='red', pointsize=40)

                ideal_left = min(ideal_left, pt[0])
                ideal_right = max(ideal_right, pt[0])
                ideal_down = min(ideal_down, pt[1])
                ideal_up = max(ideal_up, pt[1])
            else:                              # periodic solution
                pers.append(extrapolate_period(pt[0]))

                if pt[1] in RR:                # global maximum & minimum
                    ideal_down = min(ideal_down, pt[1])
                    ideal_up = max(ideal_up, pt[1])

    infls = inflection_points(fnct)
    outlist += [('infls', infls)]
    if infls not in [[], [(-oo,oo)]]:
        for pt in infls:
            if pt[0] in RR:                    # non-periodic solution
                pic += point(pt, rgbcolor='green', pointsize=20)

                ideal_left = min(ideal_left, pt[0])
                ideal_right = max(ideal_right, pt[0])
                ideal_down = min(ideal_down, pt[1])
                ideal_up = max(ideal_up, pt[1])
            else:                              # periodic solution
                pers.append(extrapolate_period(pt[0]))

                if pt[1] in RR:                # global maximum & minimum
                    ideal_down = min(ideal_down, pt[1])
                    ideal_up = max(ideal_up, pt[1])

    asmpts = asymptotes(fnct)
    outlist += [('hor_asmpts', asymptotes(fnct, astype="horizontal"))]
    outlist += [('ver_asmpts', asymptotes(fnct, astype="vertical"))]
    outlist += [('obl_asmpts', asymptotes(fnct, astype="oblique"))]
    # definitions of asymptotes
    for f in asmpts:
        # vertical asymptote - non-periodic values are pure real numbers
        if var('y') in f.args():
            ideal_left = min(ideal_left, f)
            ideal_right = max(ideal_right, f)

        # vertical asymptote - periodic values (no constant function, without x
        # in formula)
        elif var('x') not in f.args():
            pers.append(extrapolate_period(f))

        # horizontal & oblique asymptote
        else:
            ideal_down = min(ideal_down, f())
            ideal_up = max(ideal_up, f())
    ###########################################################################

    ###########################################################################
    # AUTOZOOM
    # default values for margins of graph
    l = left
    r = right
    d = down
    u = up 

    if autozoom != 0:
        # infinite range -> default range value
        if ideal_left in [-oo, oo]:
            ideal_left = left
        if ideal_right in [-oo, oo]:
            ideal_right = right
        if ideal_up in [-oo, oo]:
            ideal_up = up
        if ideal_down in [-oo, oo]:
            ideal_down = down

        hradius = abs(ideal_left - ideal_right) / 2 * autozoom
        vradius = abs(ideal_down - ideal_up) / 2

        hcenter = (ideal_left + ideal_right) / 2
        vcenter = (ideal_down + ideal_up) / 2

        if hradius != 0:
            l = hcenter - hradius
            r = hcenter + hradius
        if vradius != 0:
            d = vcenter - vradius
            u = vcenter + vradius
    ###########################################################################

    ###########################################################################
    # GENERAL FEATURES
    outlist += [('name', fnct)]

    outlist += [('even', bool(fnct == fnct(x = -x)))]
    outlist += [('odd', bool(-fnct == fnct(x = -x)))]
    outlist += [('constant', bool(self.is_real()))]     # constant function?

    undefs = undefined(fnct)
    outlist += [('undefs', undefs)]        # points where function is undefined
    for pt in undefs:
        pers.append(extrapolate_period(pt))

    the_period = oo
    for per in list(set(pers)):
        if per != 0 and fnct(x=x+per) == fnct:
            the_period = min(the_period, per)

    if the_period not in [0, oo]:
        outlist += [('period', the_period)]
        if autozoom != 0:
            l = 0
            r = the_period*autozoom
    else:
        outlist += [('period', NaN)]
    ###########################################################################

    ###########################################################################
    # RENDER PLOTS OF ASYMPTOTES & LINES & PERIODIC POINTS

    # stationary points
    if stats not in [[], [(-oo,oo)]]:
        for pt in stats:
            if pt[0] not in RR:                # periodic solution
                period = extrapolate_period(pt[0])
                if period == 0:
                    continue

                # factors (of the period) for the leftest/rightest point
                # "+1" is for addtive shift, i.e (additive_shift + z89*pi)
                leftest_factor = ceil(l/period) - 1
                rightest_factor = floor(r/period) + 1

                rng = range(leftest_factor, rightest_factor+1)
                for fctr in rng:
                    new_x = pt[0].subs( {pt[0].args()[0]:fctr} )
                    if l <= new_x <= r:
                        pic += point((new_x, fnct().subs({x:new_x})),
                                rgbcolor='red', pointsize=40)

    # tangent lines & inflections points
    if infls not in [[], [(-oo,oo)]]:
        for pt in infls:
            if pt[0] in RR:                    # non-periodic solution
                if tangents:                   # add tangents
                    tg(x) = tangent(self, pt[0])
                    if tg != []:
                        pic += plot(tg, xmin=l, xmax=r, rgbcolor='purple',
                                thickness=.5)
            else:                              # periodic solution
                period = extrapolate_period(pt[0])
                if period == 0:
                    continue

                # factors (of the period) for the leftest/rightest point
                # "+1" is for addtive shift, i.e (additive_shift + z89*pi)
                leftest_factor = ceil(l/period) - 1
                rightest_factor = floor(r/period) + 1

                rng = range(leftest_factor, rightest_factor+1)
                for fctr in rng:
                    new_x = pt[0].subs( {pt[0].args()[0]:fctr} )
                    new_y = fnct().subs({x:new_x})
                    if l <= new_x <= r and d <= new_y and new_y <= u:
                        pic += point((new_x, new_y),
                                rgbcolor='green', pointsize=20)

                        if tangents:                   # add tangents
                            tg(x) = tangent(self, new_x)
                            if tg != []:
                                pic += plot(tg, xmin=l, xmax=r,
                                        rgbcolor='purple', thickness=.5)

    # asymptotes
    for f in asmpts:
        # vertical asymptote - non-periodic values are pure real numbers
        if var('y') in f.args():
            pic += line([(f,d), (f,u)], rgbcolor='lightblue',
                    linestyle="--", thickness=2.4)

        # vertical asymptote - periodic values (no constant function, without x
        # in formula)
        elif var('x') not in f.args():
            period = extrapolate_period(pt[0])
            if period == 0:
                continue

            # factors (of the period) for the leftest/rightest point
            # "+1" is for addtive shift, i.e (additive_shift + z89*pi)
            leftest_factor = ceil(left/period) - 1
            rightest_factor = floor(right/period) + 1

            rng = range(leftest_factor, rightest_factor+1)
            for fctr in rng:
                new_x = f.subs( {f.args()[0]:fctr} )
                if left <= new_x <= right:
                    pic += line([(f,d), (f,u)], rgbcolor='lightred',
                            linestyle="--", thickness=2.4)

        # horizontal asymptote
        elif f(0) == f(1):
            pic += line([(l,f()), (r,f())], rgbcolor='lightblue',
                linestyle="--", thickness=2.4)

        # general oblique asymptote
        else:
            pic += line([(l,f().subs({f.default_variable():l})),
                (r,f().subs({f.default_variable():r}))],
                rgbcolor='lightblue', linestyle="--", thickness=2.4)
    ###########################################################################

    ###########################################################################
    # FEATURES OF 1ST DERIVATIVE
    if not self.is_real():
        if stats == []:
            end_pts = [-oo, oo]
        elif the_period not in [0, oo]:
            end_pts = []
            leftest = oo
            for pt in stats:
                if len(pt[0].args()) > 0:
                    end_pts += [pt[0].subs({pt[0].args()[0]:0})]
                    leftest = min(leftest, end_pts[-1])
            for pt in undefs:
                if len(pt.args()) > 0:
                    end_pts += [pt.subs({pt.args()[0]:0})]
                    leftest = min(leftest, end_pts[-1])
            end_pts += [leftest+the_period]
        else:                               
            end_pts = [-oo]+[pt[0] for pt in stats]+[oo]
            for pt in undefs:
                if len(pt.args()) == 0:
                    end_pts += [pt]
        end_pts = sorted(end_pts)
        outlist += [('inc', monotonicity(fnct, type="inc", hints=end_pts))]
        outlist += [('dec', monotonicity(fnct, type="dec", hints=end_pts))]
    else:
        outlist += [('inc', [(-oo,oo)])]
        outlist += [('dec', [(-oo,oo)])]
    ###########################################################################

    ###########################################################################
    # FEATURES OF 2ND DERIVATIVE
    if not self.is_real():
        if infls == [] and undefs == []:
            end_pts = [-oo, oo]
        elif the_period not in [0, oo]:
            end_pts = []
            leftest = oo
            for pt in infls:
                if len(pt[0].args()) > 0:
                    end_pts += [pt[0].subs({pt[0].args()[0]:0})]
                    leftest = min(leftest, end_pts[-1])
            for pt in undefs:
                if len(pt.args()) > 0:
                    end_pts += [pt.subs({pt.args()[0]:0})]
                    leftest = min(leftest, end_pts[-1])
            end_pts += [leftest+the_period]
        else:                               
            end_pts = [-oo]+[pt[0] for pt in infls]+[oo]
            for pt in undefs:
                if len(pt.args()) == 0:
                    end_pts += [pt]
        pts = sorted(end_pts)
        outlist += [('convex', curvature(fnct, type="convex", hints=pts))]
        outlist += [('concave', curvature(fnct, type="concave", hints=pts))]
    else:
        outlist += [('convex', [(-oo,oo)])]
        outlist += [('concave', [(-oo,oo)])]
    ############################################################################

    pl = plot(fnct, xmin=l, xmax=r, ymin=d, ymax=u, rgbcolor='black',
            detect_poles=True)
    pic += pl
    outlist += [('pic', pic)]
    forget(x, 'real')

    return dict(outlist)

def print_investigations(self, left=-2*pi, right=2*pi, down=-15, up=15,
        tangents=True, autozoom=0):
    r"""
    This function prints result of ``investigations`` to console.

    INPUT:
        
    - ``self`` - the input real-valued function ``f``
        This should be a function from ``R`` numbers to ``R``
    - ``left`` - left margin of graph
    - ``right`` - right margin of graph
    - ``down`` - down margin of graph
    - ``up`` - up margin of graph
    - ``tangents`` - bool (default: True); if True, display graph with tangent
      lines at the inflection points
    - ``autozoom`` - float (default: 1.1); if True, rescale the graph by the
      factor of ``autozoom`` in order to show as much important information
      (e. g. all reasonable stationary/inflection points, vertical asymptotes
      etc.)

    OUTPUT:

    output in console / notebook

    EXAMPLES:

    Just try running::

        sage: print_investigations(myfunc5)

    AUTHORS:

    - Karel Ha (2012-06-03): initial version
    """
    info = investigations(self, left, right, down, up, tangents, autozoom)

    print "Function", info['name']

    if info['even']: print "- even function."
    if info['odd']: print "- odd function."
    if info['constant']: print "- constant function."

    if info['undefs'] != []: print "- undefined in", info['undefs']
    if info['period'] != NaN: print "- periodic with period %s."%info['period']

    print

    print "- stationary points: ", info['stats']
    print "- increasing for: ", info['inc']
    print "- decreasing for: ", info['dec']

    print

    print "- inflection points: ", info['infls']
    print "- convex for: ", info['convex']
    print "- concave for: ", info['concave']

    print

    print "- horizontal asymptotes: ", info['hor_asmpts']
    print "- vertical asymptotes: ", info['ver_asmpts']
    print "- oblique asymptotes: ", info['obl_asmpts']

    show(info['pic'])

def tex_list(header, formulae=[]):
    r"""
    This function returns string of LaTeX's source code for ``formulae`` preceded
    by text ``header`` (as in item of ``itemize`` block)

    INPUT:
        
    - ``header`` - string; text header
    - ``formulae`` - symbolic expressions

    OUTPUT:

    string of LaTeX code

    TODO:

    - overfull TeX'x hboxes - for some reason long formulae are not
      line-wrapped

    AUTHORS:

    - Karel Ha (2012-06-10): initial version
    """
    tex_str = ""
    if formulae != []:
        tex_str += "\n\\item " + header + " $\left\{"
        for i in range(len(formulae)):
            if i != 0:          # no comma before last item
                tex_str += "; "
            tex_str += latex(formulae[i])
        tex_str += r"\right\}$"
    return tex_str

def tex_investigations(self, left=-2*pi, right=2*pi, down=-15, up=15,
        tangents=True, autozoom=0, debugmode=False, viewwith=""):
    r"""
    This function typesets the result of ``investigations`` in TeX's format.

    INPUT:
        
    - ``self`` - the input real-valued function ``f``
        This should be a function from ``R`` numbers to ``R``
    - ``left`` - left margin of graph
    - ``right`` - right margin of graph
    - ``down`` - down margin of graph
    - ``up`` - up margin of graph
    - ``tangents`` - bool (default: True); if True, display graph with tangent
      lines at the inflection points
    - ``autozoom`` - float (default: 1.1); if True, rescale the graph by the
      factor of ``autozoom`` in order to show as much important information
      (e. g. all reasonable stationary/inflection points, vertical asymptotes
      etc.)
    - ``debugmode`` - boolean (default: False); if True, show the source code
    - ``viewwith`` - boolean (default: empty string); command used for viewing
      pdf file

    OUTPUT:

    pdf file

    TODO:

    - more detailed step-by-step explanations (evenness, oddity, periodicity,
      solving "derivative's" inequalities in monotonicity & curvature parts...)
    - use of temporary files for LaTeX's source code (instead of fixed
      "/tmp/SageAlphaOut.tex")
    - option for deleting such a temporary file

    AUTHORS:

    - Karel Ha (2012-06-07): initial version
    """
    info = investigations(self, left, right, down, up, tangents, autozoom)

    tex_str = "\documentclass[10 pt]{article}"
    tex_str += "\n\usepackage[pdftex]{graphicx}"
    tex_str += "\n\usepackage{amsfonts}"
    tex_str += "\n\parindent 0pt"

    tex_str += "\n\n\\begin{document}"
    tex_str += "\n\\section*{ Function $" + latex(info['name']) + "$ \\hfil }"

    img_name = "/tmp/SageAlphaPlot.png"
    info['pic'].save(img_name)
    tex_str += "\n\n\includegraphics[scale=0.6]{" + img_name + "}"

    tex_str += "\n\n\\begin{itemize}"
    if info['even']: tex_str += "\n\\item is even function"
    if info['odd']: tex_str += "\n\\item is odd function"
    if info['constant']: tex_str += "\n\\item is constant function"
    if info['period'] != NaN:
        tex_str += "\n\\item is periodic with period " + latex(info['period'])
    tex_str += tex_list('is undefined in', formulae=info['undefs'])
    tex_str += "\n\\end{itemize}"

    tex_str += "\n\nBy analyzing the sign of the first derivative $" + \
            latex(info['name'].diff(1)) + "$:"
    tex_str += "\n\n\\begin{itemize}"
    tex_str += tex_list('stationary points are', formulae=info['stats'])
    tex_str += tex_list('is increasing for intervals', formulae=info['inc'])
    tex_str += tex_list('is decreasing for intervals', formulae=info['dec'])
    tex_str += "\n\\end{itemize}"

    tex_str += "\n\nBy analyzing the sign of the second derivative $" + \
            latex(info['name'].diff(2)) + "$:"
    tex_str += "\n\n\\begin{itemize}"
    tex_str += tex_list('inflection points are', formulae=info['infls'])
    tex_str += tex_list('is convex for intervals', formulae=info['convex'])
    tex_str += tex_list('is concave for intervals', formulae=info['concave'])
    tex_str += "\n\\end{itemize}"

    fnc = info['name']()
    if info['hor_asmpts'] != []:
        tex_str += "\n\nSince "
        tex_str += "\n\\begin{itemize}"
        l = lim(fnc, x = -oo)
        if l not in bad_limits:
            tex_str += "\n\item $\\displaystyle\lim_{x\\rightarrow-\infty}" + \
                    latex(fnc) + "=" + latex(l) + "$"
        l = lim(fnc, x = oo)
        if l not in bad_limits:
            tex_str += "\n\item and since $\\displaystyle\lim_{x\\rightarrow\infty}" + \
                    latex(fnc) + "=" + latex(l) + "$"
        tex_str += tex_list('horizontal asymptotes are', formulae=info['hor_asmpts'])
        tex_str += "\n\\end{itemize}"

    if info['ver_asmpts'] != []:
        tex_str += "\n\nSince "
        tex_str += "\n\\begin{itemize}"
        for u in info['undefs']:
            l = lim(fnc, x = u, dir='-')
            if l in [-oo, oo]:
                tex_str += "\n\item $\\displaystyle\lim_{x\\rightarrow" + \
                        latex(u) + "-}" + latex(fnc) + "=" + latex(l) + "$"
            l = lim(fnc, x = u, dir='+')
            if l in [-oo, oo]:
                tex_str += "\n\item $\\displaystyle\lim_{x\\rightarrow" + \
                        latex(u) + "+}" + latex(fnc) + "=" + latex(l) + "$"
        tex_str += tex_list('vertical asymptotes are', formulae=info['ver_asmpts'])
        tex_str += "\n\\end{itemize}"

    if info['obl_asmpts'] != []:
        tex_str += "\n\nSince "
        tex_str += "\n\\begin{itemize}"
        for inft in [-oo, oo]:
            a = lim(fnc/x, x = inft)
            if a not in bad_limits and a != 0:
                b = lim(fnc - a*x, x = oo)
                if b not in bad_limits:
                    tex_str += "\n\item $\\displaystyle\lim_{x\\rightarrow" + \
                            latex(inft) + "}\\frac{" + latex(fnc) + "}{x}"
                    tex_str += " = \\displaystyle\lim_{x\\rightarrow" + \
                            latex(inft) + "}" + latex((fnc/x).full_simplify())
                    tex_str += " = " + latex(a) + "$"
                    tex_str += "\n\item and $\\displaystyle\lim_{x\\rightarrow" + \
                            latex(inft) + "}" + latex(fnc) + "-" + latex(a) + "x"
                    tex_str += " = \\displaystyle\lim_{x\\rightarrow" + \
                            latex(inft) + "}" + latex((fnc-a*x).full_simplify())
                    tex_str += " = " + latex(b) + "$"
        tex_str += tex_list('oblique asymptotes', formulae=info['obl_asmpts'])
        tex_str += "\n\\end{itemize}"

    tex_str += "where $k \in\mathbb{Z}$."
    tex_str += "\n\\end{document}"

    # write output to TeX source code file
    fbasename = "/tmp/SageAlphaOut"
    fname = fbasename + ".tex"
    fo = open(fname, "w")
    fo.write(tex_str)
    fo.close()

    # ...compile
    os.system("pdflatex -interaction=batchmode -output-directory=/tmp " + fname)

    # ...and show it
    if viewwith != '': os.system(viewwith + ' ' + fbasename + ".pdf &")

    # show directly the TeX's source code (in console)
    if debugmode: print tex_str

def texie(self):
    r"""
    Shortcut/synonym for ``tex_investigations`` with ``okular`` as default
    viewer and ``debugmode`` turned on.

    INPUT:
        
    - ``self`` - the input real-valued function ``f``
        This should be a function from ``R`` numbers to ``R``

    OUTPUT:

    pdf file

    TODO:

    - other altenative versions (e. g. with Acrobat Reader, Evince...)

    AUTHORS:

    - Karel Ha (2012-08-15): initial version
    """
    tex_investigations(self, debugmode=True, viewwith="okular")
