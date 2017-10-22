# circle-pack

This little program is designed for visualizing the space of complex projective tori that admit a specific (square) two-circle packing. It is part of a large project to understand circle packings on complex projective structures.

## Prerequisites and running

This program requires python 3 and PyQT 5. Simply run 

`python3 circle-pack.py`

## Details

The circle packing in question corresponds to the two-circle packing on the square torus. One can see this packing on the universal cover of the square torus when the program is first started (or by hitting the `Reset` button). Each circle has its own color. On the right-hand side, the green plot shows a part of the parameter space of complex projective tori that admit this particular circle packing. The red dot is your position in the parameter space. You can drag the point or you can input parameters by hand into the fields. The `s,t` parameters are special shearing parameters that encode the complex projective torus.

The display on the left shows the developing image of the circle packing on the associated universal cover. You can change the perspective to the `exp(cz)` view, which maps one of the fixed points of the holonomy to infinity. One can also turn on dual circles.

The information on the bottom right provides the holonomy traces alongside the underlying conformal structure, given as a complex number `$\tau$` such that the torus is associated to `$\langle 1, \tau \rangle$`.

## Interface

* The red dot represents two (shearing) parameters used to compute the projective structure. You can move it! The whole parameter space is not invisible as it extends out to infinity. The space is `(s,t)` where `0 < t` and `0 < s < 2 t / (1 + 2 t^2)`. The region displayed has `0 < t < 2.5`.

* You can pan around the circle diagram on the left and use the zoom buttons in the lower right hand corner

* You can enter parameters by hand using both decimal format or (sage) syntax expressions involving `sqrt`, `exp`, `sin`, `cos` and `pi`. Your input will fail if it doesn’t land in the displayed piece of the parameter space.
    * For example : `1/(2 + sqrt(2))` is a valid expression
    * For example : `1/3^2` or `1/3**2` recall that python uses `**` for power and sage `^`. Both will work here.
    * For example : `2/(1+2*sqrt(2))` don’t forget the multiplication signs!


* You can increment and detriment the horizontal and vertical number of circles in the top right corner. If you try to decrement too far, it will reset to the original value.

* The two black curves on the parameter space correspond two families of control tori.
    * The horizontal line corresponds to rectangular tori (i.e. `$\tau$` is purely imaginary)
    * The vertical curve corresponds to tori with `$|\tau| = 1$`.

## Some pretty parameters 

* `s = 0.412, t = 0.6485`
* `s = 0.40088, t = 0.5691`
* `s = 0.580306, t = 1/sqrt(2)`
* `s = 0.45616, t = 0.58936`
* `s = 0.25, t = 1.0` in `exp(cz)` view
* `s = 0.40054, t = 1/sqrt(2)`
