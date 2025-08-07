package schoof

import (
	"goschoof/ec"
	"log"
	"math/big"
)

// ResolvePolynomialDivisionL2 naive implementation for l = 2
// for every possible x < p (the modulo of the curve) restricting to the finite field
// get all points that have y=0
// if any error happens, nil will be returned
// bad for performance, only for testing purposes
func ResolvePolynomialDivisionL2(curve *ec.EllipticCurve) []*ec.Point {
	var points []*ec.Point
	zero := big.NewInt(0)

	//for x = 0; x < p; x++
	for x := big.NewInt(0); x.Cmp(curve.GetP()) == -1; x.Add(x, big.NewInt(1)) {
		x3 := new(big.Int).Exp(x, big.NewInt(3), curve.GetP()) // x^3 mod p
		ax := new(big.Int).Mul(curve.GetA(), x)                //a * x
		ax.Mod(ax, curve.GetP())                               // a * x mod p

		val := new(big.Int).Add(x3, ax) // x^3 + a*x
		val.Add(val, curve.GetB())      //x^3 + a*x + b
		val.Mod(val, curve.GetP())      // .. mod p

		//if val == 0 (i.e. y^2 = 0 => y = 0)
		if val.Cmp(zero) == 0 {
			y := big.NewInt(0)
			p, err := ec.NewPoint(x, y)
			if err != nil {
				return nil
			}
			points = append(points, p)
		}
	}
	return points
}

// psi3 Computes ψ3(x) mod p,
// uses : 3x⁴ + 6a*x² + 12b*x - a² mod p formula
func psi3(x *big.Int, ec *ec.EllipticCurve) *big.Int {
	three := big.NewInt(3)
	six := big.NewInt(6)
	twelve := big.NewInt(12)

	x2 := new(big.Int).Exp(x, big.NewInt(2), ec.GetP())  // x^2 mod p
	x4 := new(big.Int).Exp(x2, big.NewInt(2), ec.GetP()) // x^4 mod p

	threeX4 := new(big.Int).Mul(three, x4)                               // 3x^4
	sixAx2 := new(big.Int).Mul(six, new(big.Int).Mul(ec.GetA(), x2))     // 6*ax^2
	twelveBx := new(big.Int).Mul(twelve, new(big.Int).Mul(ec.GetB(), x)) // 12*b*x
	a2 := new(big.Int).Exp(ec.GetA(), big.NewInt(2), ec.GetP())          // a^2 mod p

	val := new(big.Int).Add(threeX4, sixAx2) // 3x^4 + 6*ax^2
	val.Add(val, twelveBx)                   // .. + 12*b*x
	val.Sub(val, a2)                         // .. - a^2
	val.Mod(val, ec.GetP())                  // .. mod p
	return val

}

// ResolvePolynomialDivisionL3 naive implementation for l = 3
// gets all the x's that fulfills: ψ3(x) mod p ≡ 0 mod p
// then try to get all the y for all the x's by resolving Weierstrass equation for the given x
// bad for performance, only for testing purposes
func ResolvePolynomialDivisionL3(curve *ec.EllipticCurve) []*ec.Point {
	var points []*ec.Point

	var xs []*big.Int
	//for x = 0; x < p; x++
	for x := big.NewInt(0); x.Cmp(curve.GetP()) == -1; x.Add(x, big.NewInt(1)) {
		psi3 := psi3(x, curve)
		//if ψ3(x) mod p ≡ 0 mod p
		if psi3.Cmp(big.NewInt(0)) == 0 {
			xp := new(big.Int).Set(x)
			xs = append(xs, xp)
		}
	}

	//compute y coordinates for all x's, appends to list
	for _, x := range xs {
		y, ok := curve.ProcessYFrom(x)
		if !ok {
			continue
		}
		point, err := ec.NewPoint(x, y)
		if err != nil {
			log.Fatalf(err.Error())
		}
		points = append(points, point)
	}

	return points
}

// PSI_l
// if l is odd $ ψ_{2m+1} = ψ_{m+2} * ψ³_{m} − ψ_{m−1] * ψ³_{m+1} $
// else $ψ_{2m}=\frac{ψ_m}{2y}(ψ_{m+2} * ψ²_{m−1}−ψ_{m−2} * ψ²_{m+1})$
func PSI_l(curve *ec.EllipticCurve, l int64) []*ec.Point {
	if l == 2 {
		return ResolvePolynomialDivisionL2(curve)
	}
	if l == 3 {
		return ResolvePolynomialDivisionL3(curve)
	}
	//TODO: add a cache system to prevent repetitive calculations
	return nil
}
