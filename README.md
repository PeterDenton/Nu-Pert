Nu-Pert
=

[![Build Status](https://travis-ci.org/PeterDenton/Nu-Pert.svg?branch=master)](https://travis-ci.org/PeterDenton/Nu-Pert)

## Overview
This code provides a very precise, yet still useful analytically, means of calculating the neurtino transition probabilities in uniform matter densities.
We have included our results to various orders in a perturbative expansion.
Our perturbative expansion is proprotional to the ratio of the two measured <nobr>Delta m<sup>2</sup>'s</nobr>.
It is smaller than 1.5% and is exactly zero in vacuum.

Zeroth order is more than sufficient for any modern neutrino experiment with precision on the order of 10<sup>-4</sup>.
First order improves this by about three more orders of magnitude.
Second order is also included through a different code that provides about three more orders of magnitude in precision, and, in some cases, hits the limit of double precision.
These results contain all orders in the matter potential and in s<sub>13</sub>.

## About the code
Equation numbers are provided where possible which reference [arXiv:1604.xxxxx](https://arxiv.org/abs/1604.xxxxx) unless otherwise stated.
The various namespaces throughout the code refer to the various steps and means of calculating the transition probabilities.

First, the `Hat` namespace is in reference to the intermediate basis after the 23 and 13 rotations.
The `Hat` namespace refers to the results from the [arXiv:1505.01826](https://arxiv.org/abs/1505.01826) paper by S. Parke and H. Minakata.
This namespace contains the eigenvalues and the angle phi.
This namespace also includes the explicit formulas for P(nu<sub>e</sub>-&gt;nu<sub>e</sub>) and P(nu<sub>e</sub>-&gt;nu<sub>mu</sub>) in this basis.

Next is the `Check` namespace which is in reference to the final basis after an additional 12 rotation.
This namespace contains the final eigenvalues at zeroth and second order (zeroth order also contains all first order corrections), and the angle psi.
It also contains the W and V matrices as well as a transition probability calculation.
The probability calculated in this namespace is based on multiplying out the V matrix.
This is identical to the `GF` namespace method at zeroth order, but at first order it contains some higher order terms.
This method is coded through second order, but could be expanded to higher orders.

The `GF` namespace refers to the general form of the first order transition probabilities presented in eqs. blah, and explicitly shown in table 1.
This namespace provides an order-by-order calculation for zeroth or first order in that the first order calculation provides no manifestly higher order terms.

The vacuum values of the oscillation parameters are described in Constant.cpp.
Since there is no consistent way for all channels to switch the mass ordering, we take the simple choice of changing the sign on Delta m<sup>2</sup><sub>31</sub>.

See Examples.cpp for several examples.

## Some further notes
The antineutrino transition probabilities are the same as the neutrino transition probabilities except with E-&gt;-E.
That is, L/E-&gt;-L/E and a-&gt;-a.

P(nu<sub>a</sub>-&gt;nu<sub>b</sub>) is related to P(nu<sub>b</sub>-&gt;nu<sub>a</sub>) by sending L-&gt;-L.

## Support
If you have questions or encounter any problems when running *Nu-Pert*, please use github's [issue tracker](https://github.com/PeterDenton/Nu-Pert/issues).

This code is free to use, copy, distribute, and modify.
If you use this code or any modification of this code, we request that you reference the relevant publication, [arXiv:1604.xxxxx](https://arxiv.org/abs/1604.xxxxx).
