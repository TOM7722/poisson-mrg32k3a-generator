# Poisson Distribution Generator using Inversion Method

## Detailed Description
This project implements an inversion method to generate random variables following a Poisson distribution with parameter λ. The uniform numbers generation is performed using the MRG32K3a generator.

MRG32K3a is a combined pseudo-random number generator developed by Pierre L'Ecuyer (UdeM).

### How MRG32K3a Works
First component: xn = (1403580 * x[n-2] - 810728 * x[n-3]) mod m1
where m1 = 4294967087
Second component: yn = (527612 * y[n-1] - 1370589 * y[n-3]) mod m2
where m2 = 4294944443

Final combination: un = ((xn - yn) mod m1) / m1

## Inversion Method
The implemented method uses the principle of inverting the cumulative distribution function:
1. Generation of a uniform number U on [0,1] using MRG32K3a
2. Search for the smallest value k such that F(k) ≥ U
   where F is the cumulative distribution function of the Poisson distribution

Two approaches are proposed:
- Sequential search (simpler but less efficient)
- Indexed search (faster for large simulations)

### Required Libraries
```python
import math            # Built-in mathematical functions
import matplotlib.pyplot as plt   # Version 3.7.1 or higher
import numpy as np              # Version 1.21.0 or higher
from scipy.stats import poisson  # Version 1.7.0 or higher
'''

### Source
L'Ecuyer, P. (1999). "Good Parameters and Implementations for Combined Multiple Recursive Random Number Generators". Operations Research, 47(1), 159-164.
[DOI: 10.1287/opre.47.1.159](https://doi.org/10.1287/opre.47.1.159)
