package schoof

import (
	"goschoof/ec"
	"goschoof/polynom"
	"goschoof/utils"
	"log"
	"math/big"
)

// PSICache Cache for PSI calculations
// (storing result and preventing deep recursion for already done calculations).
type PSICache struct {
	curve *ec.EllipticCurve
	cache map[string]*polynom.Polynom
}

func NewPSICache(curve *ec.EllipticCurve) *PSICache {
	return &PSICache{
		curve: curve,
		cache: make(map[string]*polynom.Polynom),
	}
}

func CountTorsion2PointsFromPoly(curve *ec.EllipticCurve) int {
	psi2 := BuildPolynomL2(curve)

	p := curve.GetP()
	if p.BitLen() > 32 { // if p is big
		return 1
	}

	count := 0
	zero := big.NewInt(0)

	for x := big.NewInt(0); x.Cmp(p) < 0; x.Add(x, big.NewInt(1)) {
		result := evaluatePolyAt(psi2, x)
		if result.Cmp(zero) == 0 {
			count++
		}
	}

	return count
}

// evaluatePolyAt
func evaluatePolyAt(poly *polynom.Polynom, x *big.Int) *big.Int {
	if len(poly.Coefficients) == 0 {
		return big.NewInt(0)
	}

	result := big.NewInt(0)
	xPower := big.NewInt(1) // x^0 = 1

	for i, coeff := range poly.Coefficients {
		if i > 0 {
			xPower.Mul(xPower, x)
			xPower.Mod(xPower, poly.P)
		}

		term := new(big.Int).Mul(coeff, xPower)
		term.Mod(term, poly.P)

		result.Add(result, term)
		result.Mod(result, poly.P)
	}

	return result
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
	// for x = 0; x < p; x++
	for x := big.NewInt(0); x.Cmp(curve.GetP()) == -1; x.Add(x, big.NewInt(1)) {
		psi3 := psi3(x, curve)
		// if ψ3(x) mod p ≡ 0 mod p
		if psi3.Cmp(big.NewInt(0)) == 0 {
			xp := new(big.Int).Set(x)
			xs = append(xs, xp)
		}
	}

	// compute y coordinates for all x's, appends to list
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

func BuildPolynomL4(curve *ec.EllipticCurve) *polynom.Polynom {
	// ψ₄ = 2(4x⁶ + 5ax⁴ + 20bx³ - 5a²x² - 4abx - 8b² - a³)

	a := curve.GetA()
	b := curve.GetB()
	p := curve.GetP()

	a2 := new(big.Int).Exp(a, big.NewInt(2), p) // a²
	a3 := new(big.Int).Exp(a, big.NewInt(3), p) // a³
	b2 := new(big.Int).Exp(b, big.NewInt(2), p) // b²
	ab := new(big.Int).Mul(a, b)                // ab
	ab.Mod(ab, p)

	coeff0 := new(big.Int).Mul(big.NewInt(-8), b2)
	coeff0.Sub(coeff0, a3)
	coeff0.Mul(coeff0, big.NewInt(2))
	coeff0.Mod(coeff0, p)

	coeff1 := new(big.Int).Mul(big.NewInt(-8), ab)
	coeff1.Mod(coeff1, p)

	coeff2 := new(big.Int).Mul(big.NewInt(-10), a2)
	coeff2.Mod(coeff2, p)

	coeff3 := new(big.Int).Mul(big.NewInt(40), b)
	coeff3.Mod(coeff3, p)

	coeff4 := new(big.Int).Mul(big.NewInt(10), a)
	coeff4.Mod(coeff4, p)

	coeff5 := big.NewInt(0)

	coeff6 := big.NewInt(8)

	poly := polynom.NewPolynom([]*big.Int{
		coeff0, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6,
	}, curve.GetP())

	return poly
}

// PSI_l - computes psi_l
func PSI_l(curve *ec.EllipticCurve, l *big.Int, psiCache *PSICache) *polynom.Polynom {
	if l.Sign() <= 0 {
		log.Panicf("psi index (%s) <= 0: forbidden", l.String())
	}

	lStr := l.String()
	// short-circuiting existing calculated psi values
	if poly, ok := psiCache.cache[lStr]; ok {
		return poly
	}

	var res *polynom.Polynom

	switch {
	case l.Cmp(big.NewInt(0)) == 0:
		res = polynom.NewPolynom([]*big.Int{big.NewInt(0)}, curve.GetP())
		psiCache.cache[lStr] = res
		return res
	case l.Cmp(big.NewInt(1)) == 0:
		// if l = 1, the corresponding polynomial is only 1, with no x
		res = polynom.NewPolynom([]*big.Int{big.NewInt(1)}, curve.GetP())
		psiCache.cache[lStr] = res
		return res
	case l.Cmp(big.NewInt(2)) == 0:
		res = BuildPolynomL2(curve)
		psiCache.cache[lStr] = res
		return res
	case l.Cmp(big.NewInt(3)) == 0:
		res = BuildPolynomL3(curve)
		psiCache.cache[lStr] = res
		return res
	case l.Cmp(big.NewInt(4)) == 0:
		res = BuildPolynomL4(curve)
		psiCache.cache[lStr] = res
		return res
	}

	//l is odd: $ ψ_{2m+1} = ψ_{m+2} * ψ³_{m} − ψ_{m−1] * ψ³_{m+1} $
	// so
	// 	   2m+1 = l
	//<==> 2m = l - 1
	//<==> m = (l - 1) / 2
	if new(big.Int).Mod(l, big.NewInt(2)).Cmp(big.NewInt(0)) != 0 {
		lMinus1 := new(big.Int).Sub(l, big.NewInt(1))
		m := new(big.Int).Div(lMinus1, big.NewInt(2))
		log.Printf("psi_l:: l=%v, current m=%v :: lMinus1=%v (not launching recursion yet)", l, m, lMinus1)

		// to have ψ_{2m+1}, we have to calculate ψ_{m+2}, ψ³_{m}, ψ_{m-1} & ψ³_{m+1}
		// m will be recursively a 'l' value
		mPlus2 := new(big.Int).Add(m, big.NewInt(2))
		log.Printf("psi_l:: l=%v, current m=%v :: mPlus2=%v (about to start recursion for PsiMp2, with mPlus2(%v) as l)", l, m, mPlus2, mPlus2)
		PsiMp2 := PSI_l(curve, mPlus2, psiCache)
		log.Printf("psi_l:: l=%v, current m=%v :: PsiMp2=%v (finished recursion for PsiMp2)", l, m, PsiMp2)

		log.Printf("psi_l:: l=%v, current m=%v :: (about to start recursion for PsiM, with m(%v) as l)", l, m, m)
		PsiM := PSI_l(curve, m, psiCache)
		log.Printf("psi_l:: l=%v, current m=%v :: PsiM=%v (finished recursion for PsiM)", l, m, PsiM)
		PsiMcubed := PsiM.Mul(PsiM).Mul(PsiM)

		mMinus1 := new(big.Int).Sub(m, big.NewInt(1))
		PsiMm1 := PSI_l(curve, mMinus1, psiCache)

		mPlus1 := new(big.Int).Add(m, big.NewInt(1))
		PsiMp1 := PSI_l(curve, mPlus1, psiCache)
		PsiMp1cubed := PsiMp1.Mul(PsiMp1).Mul(PsiMp1)

		term1 := PsiMp2.Mul(PsiMcubed)   // ψ_{m+2} * ψ³_{m}
		term2 := PsiMm1.Mul(PsiMp1cubed) // ψ_{m−1] * ψ³_{m+1}
		res = term1.Sub(term2)           // ψ_{m+2} * ψ³_{m} − ψ_{m−1] * ψ³_{m+1}
	} else {
		// forbidden even l
		log.Panicf("psi_l with an even number is not allowed: %v", l)
	}

	psiCache.cache[lStr] = res
	return res
}

func getSmallL(curve *ec.EllipticCurve) []*big.Int {
	sqrtp := new(big.Int).Sqrt(curve.GetP())         // floor(sqrt(p))
	target := new(big.Int).Mul(big.NewInt(4), sqrtp) // 4*sqrt(p)

	M := big.NewInt(1)
	var ls []*big.Int

	for l := int64(3); ; l += 2 { // no even numbers
		if !utils.IsPrime(l) {
			continue
		}
		lBig := big.NewInt(l)
		// i l == p skip
		if lBig.Cmp(curve.GetP()) == 0 {
			continue
		}

		tmp := new(big.Int).Mul(M, lBig)

		if tmp.Cmp(target) > 0 {
			break
		}
		M.Set(tmp)
		ls = append(ls, new(big.Int).Set(lBig))
	}
	return ls
}

func crtUpdate(T, M, c, l *big.Int) (*big.Int, *big.Int) {
	// k = ((c - T) mod l) * (M^{-1} mod l) mod l
	diff := new(big.Int).Sub(c, new(big.Int).Mod(T, l))
	diff.Mod(diff, l)
	if diff.Sign() < 0 {
		diff.Add(diff, l)
	}

	Minv := new(big.Int).ModInverse(new(big.Int).Mod(M, l), l)
	if Minv == nil {
		panic("gcd(M,l) != 1")
	}

	k := new(big.Int).Mul(diff, Minv)
	k.Mod(k, l)

	Tnew := new(big.Int).Add(T, new(big.Int).Mul(k, M))
	Mnew := new(big.Int).Mul(M, l)
	return Tnew, Mnew
}

func Schoof(curve *ec.EllipticCurve) *big.Int {
	ls := getSmallL(curve)
	T := big.NewInt(0) // t mod M
	M := big.NewInt(1) // prod of ℓ treated
	target := new(big.Int).Mul(big.NewInt(4), new(big.Int).Sqrt(curve.GetP()))

	// l=2 ψ₂
	count2 := CountTorsion2PointsFromPoly(curve)
	t2 := big.NewInt(int64(count2 % 2)) // t mod 2
	T, M = crtUpdate(T, M, t2, big.NewInt(2))

	x := polynom.NewPolynom([]*big.Int{big.NewInt(0), big.NewInt(1)}, curve.GetP())

	cache := NewPSICache(curve)
	for _, l := range ls {
		log.Printf("schoof::Schoof > treating l=%d", l)
		h := PSI_l(curve, l, cache)
		Xp := x.PowMod(curve.GetP(), h)
		Xp2 := Xp.PowMod(curve.GetP(), h)
		c := computeTmodL(curve, l, h, x, Xp, Xp2, curve.GetP(), cache)
		T, M = crtUpdate(T, M, c, l)
		if M.Cmp(target) > 0 {
			break
		}
	}

	N := new(big.Int).Add(curve.GetP(), big.NewInt(1))
	N.Sub(N, T)
	return N
}

func computeTmodL(curve *ec.EllipticCurve, l *big.Int, h *polynom.Polynom, xPoly, Xp, Xp2 *polynom.Polynom, p *big.Int, cache *PSICache) *big.Int {
	// loop on c = 0..l-1
	zero := big.NewInt(0)
	one := big.NewInt(1)
	for c := new(big.Int).Set(zero); c.Cmp(l) < 0; c.Add(c, one) {
		// target = π^2 - c*π + p on x coordinate
		// π(x,y) = (x^p, y^p) so π(x) = x^p = Xp
		// π^2(x) = (x^p)^p = x^(p^2) = Xp2
		// => target = Xp2 - c*Xp + p*1

		// compute c*Xp
		comb := Xp.Scale(c)
		_, comb = comb.DivMod(h)

		target := Xp2.Sub(comb)
		pCoeff := new(big.Int).Mod(p, curve.GetP())
		if target.Coeff(0).Add(target.Coeff(0), pCoeff); true {
			target.SetCoeff(0, new(big.Int).Add(target.Coeff(0), pCoeff))
		}
		_, target = target.DivMod(h)

		if target.IsZero() {
			return new(big.Int).Set(c)
		}

		d := polynom.GCDPolynom(h, target)
		if d != nil && d.Degree() > 0 && d.Degree() < h.Degree() {
			if d.Degree() <= h.Degree()-d.Degree() {
				h = d
			} else {
				h = polynom.DivExact(h, d)
			}
			// recompute Xp, Xp2 with new h
			Xp = xPoly.PowMod(p, h)
			Xp2 = Xp.PowMod(p, h)
			c.Sub(c, one)
			continue
		}
	}
	panic("no c found!")
}
