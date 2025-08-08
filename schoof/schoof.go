package schoof

import (
	"goschoof/ec"
	"goschoof/polynom"
	"log"
	"math/big"
)

// PSICache Cache for PSI calculations
// (storing result and preventing deep recursion for already done calculations).
type PSICache struct {
	curve *ec.EllipticCurve
	cache map[int64]*polynom.Polynom
}

func NewPSICache(curve *ec.EllipticCurve) *PSICache {
	return &PSICache{
		curve: curve,
		cache: make(map[int64]*polynom.Polynom),
	}
}

// GetPointsOfL2 naive implementation for l = 2
// for every possible x < p (the modulo of the curve) restricting to the finite field
// get all points that have y=0
// if any error happens, nil will be returned
// bad for performance, only for testing purposes
func GetPointsOfL2(curve *ec.EllipticCurve) []*ec.Point {
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

func BuildPolynomL2(curve *ec.EllipticCurve) *polynom.Polynom {
	fourB := new(big.Int).Mul(big.NewInt(4), curve.GetB())
	fourA := new(big.Int).Mul(big.NewInt(4), curve.GetA())
	poly := polynom.NewPolynom([]*big.Int{
		fourB, fourA, big.NewInt(0), big.NewInt(4),
	}, curve.GetP())
	return poly
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

// BuildPolynomL3 3x⁴ + 6a*x² + 12b*x - a² mod p
func BuildPolynomL3(curve *ec.EllipticCurve) *polynom.Polynom {
	negASquared := new(big.Int).Exp(curve.GetA(), big.NewInt(2), curve.GetP()) // a²
	negASquared.Neg(negASquared)                                               // - a²

	twelveB := new(big.Int).Mul(big.NewInt(12), curve.GetB()) // 12b
	sixA := new(big.Int).Mul(big.NewInt(6), curve.GetA())     // 6a
	poly := polynom.NewPolynom([]*big.Int{
		negASquared, twelveB, sixA, big.NewInt(0), big.NewInt(3),
	}, curve.GetP())
	return poly
}

// PSI_l
// if l is odd $ ψ_{2m+1} = ψ_{m+2} * ψ³_{m} − ψ_{m−1] * ψ³_{m+1} $
// else $ψ_{2m}=\frac{ψ_m}{2y}(ψ_{m+2} * ψ²_{m−1}−ψ_{m−2} * ψ²_{m+1})$ (but will change, no division & no y)
func PSI_l(curve *ec.EllipticCurve, l int64, psiCache *PSICache) *polynom.Polynom {
	//short-circuiting existing calculated psi values
	if poly, ok := psiCache.cache[l]; ok {
		return poly
	}

	var res *polynom.Polynom

	switch l {
	case 0:
		res = polynom.NewPolynom([]*big.Int{big.NewInt(0)}, curve.GetP())
		psiCache.cache[l] = res
		return res
	case 1:
		//if l = 1, the corresponding polynomial is only 1, with no x
		res = polynom.NewPolynom([]*big.Int{big.NewInt(1)}, curve.GetP())
		psiCache.cache[l] = res
		return res
	case 2:
		res = BuildPolynomL2(curve)
		psiCache.cache[l] = res
		return res

	case 3:
		res = BuildPolynomL3(curve)
		psiCache.cache[l] = res
		return res
	case 4:
		// avoiding infinite looping with m+2 case (2+2 = 4 since m = l / 2 & for l = 4
		// would be calculating m = 4 (we need m+2) to get l = 4
		// leading to stack overflow for infinite recursion
		p2 := PSI_l(curve, 2, psiCache)
		p3 := PSI_l(curve, 3, psiCache)
		res := p2.Mul(p2).Mul(p3.Mul(p3))
		psiCache.cache[l] = res
		return res
	}

	//l is odd: $ ψ_{2m+1} = ψ_{m+2} * ψ³_{m} − ψ_{m−1] * ψ³_{m+1} $
	// so
	// 	   2m+1 = l
	//<==> 2m = l - 1
	//<==> m = (l - 1) / 2
	if l%2 != 0 {
		m := (l - 1) / 2
		// to have ψ_{2m+1}, we have to calculate ψ_{m+2}, ψ³_{m}, ψ_{m-1} & ψ³_{m+1}
		// m will be recursively a 'l' value
		PsiMp2 := PSI_l(curve, m+2, psiCache)

		PsiM := PSI_l(curve, m, psiCache)
		PsiMcubed := PsiM.Mul(PsiM).Mul(PsiM)

		PsiMm1 := PSI_l(curve, m-1, psiCache)

		PsiMp1 := PSI_l(curve, m+1, psiCache)
		PsiMp1cubed := PsiMp1.Mul(PsiMp1).Mul(PsiMp1)

		term1 := PsiMp2.Mul(PsiMcubed)   // ψ_{m+2} * ψ³_{m}
		term2 := PsiMm1.Mul(PsiMp1cubed) //ψ_{m−1] * ψ³_{m+1}
		res = term1.Sub(term2)           // ψ_{m+2} * ψ³_{m} − ψ_{m−1] * ψ³_{m+1}
	} else {
		// l is even
		// $ψ_{2m}=\frac{ψ_m}{2y}(ψ_{m+2} * ψ²_{m−1}−ψ_{m−2} * ψ²_{m+1})$
		// 	    l = 2m
		// <==> l/2 = m

		// since we cannot do ψ_m/2y (we want polynomials with only x)
		// we calculate the T_l² polynomial (not ψ_l nor ψ_l²)
		// to allow eliminating the y and the division
		m := l / 2
		M := PSI_l(curve, m, psiCache)
		A := PSI_l(curve, m+2, psiCache) // m must NOT be 2 in case the cache doesn't contain l=4, or will SO
		B := PSI_l(curve, m-1, psiCache)
		C := PSI_l(curve, m-2, psiCache)
		D := PSI_l(curve, m+1, psiCache)

		B2 := B.Mul(B)                // B²
		D2 := D.Mul(D)                // D²
		M2 := M.Mul(M)                // M²
		E := A.Mul(B2).Sub(C.Mul(D2)) // A*B² - C*D²
		res = M2.Mul(E.Mul(E))        // M² * E²
	}

	psiCache.cache[l] = res
	return res
}
