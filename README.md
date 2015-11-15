# Polynomial library

A quick & dirty C++ library to handle polynoms.
This is not intended for serious usage, this is just an example / exercise.

Poly is a template class allowing to represent polynoms with various types of coefficient. 
The test include examples with integers, boost arbitrary precision integers, rationals and finite fields (GF(n)).
The examples on finite fields are based on another class provided implementing GF(n) element for n prime (i.e simple n-modular arithmetic).

Two implementations of polynomials over GF(2) are also included providing better performance for this specific case.

The code includes tentative implementation of polynomial factorization algorithms:
 - Cantor-Zassenhaus (I'm not 100% sure it's coded properly as the probability of a random polynom to 'work' seems much lower than the theory)
 - McEliece, based on this paper: 
 - Knuth, wrongly named but based on exercise 30 of paragraph 4.6.2 of volume 2 of The Art of Computer Programming. (Doesn't seem to work as exepected either)
