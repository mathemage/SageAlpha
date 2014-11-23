SageAlpha
=========

- university project for course Individual Software Project
- 2nd year of bachelor studies

![alt text][SageMath-logo]

About
-----

Enhanced solving of easier exercises from mathematical analysis (aimed at the difficulty level of freshmen students at the Faculty of Mathematics and Physics of the Charles University in Prague, Czech Republic).

Written in Python as a SageMath's module.

Features
--------

- `investigations` descriptive investigation of a real-valued function
  - general properties (odd, even, periodic, constant)
  - properties of 1st derivative (stationary points, increasing/decreasing intervals)
  - properties of 2nd derivative (inflection points, convex/concave intervals)
  - asymptotes (horizontal, vertical & oblique, with explained calculation)
  - rendering the plot (with tangents at inflection points & autozoom option)
  - two output formats - text-based (console/notebook) & LaTeX-based (pdf document example `src/example-SageAlpha.pdf`)

Versions
--------

- *v1.0.0*  (Sep 2012)
  - tested on Sage 5.1, Sage 5.0 and Sage 4.8
  - mostly working (known exceptions: `monotonicity(myfunc7)` ; `investigations(x^5-x^3-x^2+x+1)` )

TODO
----

- `undefined` further causes of *undefinition*: negative arguments of square roots (and their fourth, sixth root analogies...), *arcsin* & *arccos* arguments...
- `stationary_points` points for x-values where 1st derivative is not defined
- `inflection_points` dtto
- `asymptotes` verify proper functionality of periodical (vertical) asymptotes - e. g. in *tan(x)*
- `monotonicity`
  - **either** improve splitting points by adding `stationary_points` where 1st derivative is undefined
  - **or** find better way to locate these intervals (preferably using `solve` that will find intervals of positive/negative 1st derivative)
- `curvature` dtto
- `investigations`
  - simplify complex formulae of points for rendering in the plot (see TODO example in the docstring via `investigations?`)
  - beta-tests of periodical (vertical) asymptotes' rendering
  - more comprehensible comments/documentation for *autozoom* feature
- `tex_list` line-wrapping of long formulae (or when there are simply too many of them, as well)
- `tex_investigations`
  - more detailed step-by-step explanations (see TODO part in the docstring via `tex_investigations?`)
  - use of temporary files for LaTeX's source code (instead of fixed `/tmp/SageAlphaOut.tex`)
  - option for deleting such a temporary file
- `texie` other altenative versions (e. g. with Acrobat Reader, Evince...)
- further functionalities - e. g.:
  - step-by-step solving of simpler integrals
  - convergence/divergence of easier power series (with tutorials on which criteria to use)
  - multidimensional stuff (implicit functions, Hessians, Lagrangians, ...)
  - ...


[SageMath-logo]: ./sage-squared.png "SageAlpha is a module in SageMath"
