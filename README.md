# goschoof : Schoof algorithm implementation in Golang

### project still in progress, thus, the Schoof implementation is not done yet, & will be followed by a Schoof-Elkies-Atkin implementation when done

This implementation also provides an elliptic curves implementation (probably not the best optimized) focused on cryptography 
(by checking for groups, primality of $p$ for the elliptic curve, non-singularity, etc.).


## Elliptic curves

> Defined with the following equation (Weierstrass equation):
> 
> $ y^2 = x^3 + Ax + B$ with $A,B\in\mathbb F_q$, $q = p^n$, $p$ a prime number and $n$ an integer $\ge 1$, $p\neq2,3$.

## Polynomials
To get points of a given order (ℓ) of the curve.

### for $ℓ = 2$
All points $(x,y)$ of the curve having $y=0$, $\forall x < p$.

### for $ℓ = 3$
Test all x with this formula:

$3x⁴ + 6a x² + 12bx - a² \mod p \iff \psi_3$ 

We consider $\psi_x$ a polynomial, with $x$ a . 

Then, retrieve all points of the curve (i.e. get all the y's that make a $(x,y)$ pair of 
coordinates that are points of the curve). 

### for $ℓ \geq 4$ (general case)
Recursive division polynomials generation

If $ℓ \mod 2 \neq0$

$ψ_{2m+1} = ψ_{m+2} \times ψ³_{m} − ψ_{m−1} \times ψ³_{m+1}$

If $ℓ \mod 2 =0$

$ψ_{2m}=\frac{ψ_m}{2y}(ψ_{m+2} \times ψ²_{m−1}−ψ_{m−2} \times ψ²_{m+1})$

## References

- Hasse theorem [Wikipedia](https://en.wikipedia.org/wiki/Hasse%27s_theorem_on_elliptic_curves)

> $$|N - (q+1)| \le 2\sqrt{q}$$
> 
> with $N$ the number of points of the elliptic curve, on a finite field with $q$ elements. 
- Schoof-Elkies-Atkin algorithm (extension of Schoof algorithm) https://en.wikipedia.org/wiki/Schoof%E2%80%93Elkies%E2%80%93Atkin_algorithm
- Schoof's algorithm for Counting Points on $E(\mathbb F_q)$, Gregg Musiker, 7th December 2005 https://www-users.cse.umn.edu/~musiker/schoof.pdf
- 
- `secp256k1` elliptic curve (used for bitcoin) https://en.bitcoin.it/wiki/Secp256k1

<img src="assets/Secp256k1.png" style="background: aliceblue">