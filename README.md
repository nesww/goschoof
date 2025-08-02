# goschoof : Schoof algorithm implementation in Golang

This implementation also provides an elliptic curves implementation (probably not the best optimized) focused on cryptography 
(by checking for groups, primality of $p$ for the elliptic curve, non-singularity, etc.).


## Elliptic curves

> Defined with the following equation (Weierstrass equation):
> 
> $$ y^2 = x^3 + Ax + B$$ with $A,B\in\mathbb F_q$, $q = p^n$, $p$ a prime number and $n$ an integer $\ge 1$, $p\neq2,3$.

## References

- Hasse theorem [Wikipedia](https://en.wikipedia.org/wiki/Hasse%27s_theorem_on_elliptic_curves)

> $$|N - (q+1)| \le 2\sqrt{q}$$
> 
> with $N$ the number of points of the elliptic curve, on a finite field with $q$ elements. 
- Using Miller Rabin test for prime numbers implementation from : [miller-rabin](https://github.com/boppreh/miller-rabin), slightly modified to fit the usage and module import.
- `secp256k1` elliptic curve (used for bitcoin) https://en.bitcoin.it/wiki/Secp256k1

<img src="assets/Secp256k1.png" style="background: aliceblue">