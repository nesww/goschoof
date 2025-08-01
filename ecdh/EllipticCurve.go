package ecdh

import (
	"fmt"
	"math/big"
)

type EllipticCurve struct{
	a *big.Int
	b *big.Int
	p *big.Int //has to be prime
}

func NewEllipticCurve(a, b, p *big.Int) (EllipticCurve, error) {
	if !MillerRabin(p) {
		return EllipticCurve{}, fmt.Errorf("Given p %s for prime number as modulo of the curve is not prime", p)
	}

	modA := new(big.Int).Mod(a,p)
	modB := new(big.Int).Mod(b,p)
	modP := new(big.Int).Set(p)

	return EllipticCurve{modA,modB,modP}, nil;
}

// True iff the elliptic curve isn't singular i.e. not crossing itself and not pointed, 
// reducing risks of cryptographic vulnerability.
func (ec *EllipticCurve) IsNonSingular() bool {
	// 4 * a ** 3
	four := big.NewInt(4);
	aCubed := new(big.Int).Exp(ec.a, big.NewInt(3), nil)
	q := new(big.Int).Mul(four, aCubed);

	// 27 * b ** 2
	twentySeven := big.NewInt(27);
	bSquared := new(big.Int).Exp(ec.b, big.NewInt(2), nil)
	d := new(big.Int).Mul(twentySeven, bSquared)

	q.Add(q,d); // 4 * a**3  +  27 * b**2

	q.Mod(q, ec.p) // q = q % p

	return q.Cmp(big.NewInt(0)) != 0
}

// Calculates the slope existing between the segment linking both given points.
// Takes into consideration that P could be Q.
func (ec *EllipticCurve) CalcSlope(p *Point, q *Point) *big.Int {
	if p.Equals(q) {
		// (3 * (p.x^2)  + a)
		three := big.NewInt(3); // 3
		pxSquared := new(big.Int).Exp(p.x, big.NewInt(2), nil) // p.x^2
		left := new(big.Int).Mul(three, pxSquared) // 3 * p.x^2
		left.Add(left, ec.a); // 3 * p.x^2 + a
		left.Mod(left, ec.p); // .. mod p

		// modular_inverse of 2 * p.y % p
		two_py := new(big.Int).Mul(big.NewInt(2), p.y) // 2 * p.y
		right := new(big.Int).ModInverse(two_py, ec.p) // mod_inverse of 2 * p.y mod p

		//in case of no modular_inverse (rare)
		if right == nil {
			return nil;
		}

		left.Mul(left, right); // (3 * (p.x)^2 + a) * ( mod_inv(2 * p.y) mod p)
		left.Mod(left, ec.p) // ^^^^^^ mod p
		return left;


	} else {
		// P != Q

		//left
		qy_Sub_py := new(big.Int).Sub(q.y, p.y); // q.y - p.y

		//right
		// modular_inverse of q.x - p.x mod ec.p
		qx_Sub_px := new(big.Int).Sub(q.x, p.x);
		right := new(big.Int).ModInverse(qx_Sub_px, ec.p); // mod_inverse of (q.x - p.x) mod ec.p
		
		//if no modular inverse (rare)
		if right == nil {
			return nil;
		}
		
		qy_Sub_py.Mul(qy_Sub_py, right); // (q.y - p.y) * mod_inverse of (q.x - p.x) mod ec.p
		qy_Sub_py.Mod(qy_Sub_py, ec.p) // ^^^^^^ mod ec.p

		return qy_Sub_py;
	}
}

