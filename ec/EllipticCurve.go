package ec

import (
	"fmt"
	"goschoof/utils"
	"log"
	"math/big"
)

type EllipticCurve struct {
	a *big.Int
	b *big.Int
	p *big.Int //has to be prime
}

func (ec *EllipticCurve) GetA() *big.Int {
	return ec.a
}

func (ec *EllipticCurve) GetB() *big.Int {
	return ec.b
}

func (ec *EllipticCurve) GetP() *big.Int {
	return ec.p
}

func NewEllipticCurve(a, b, p *big.Int) (*EllipticCurve, error) {
	if !utils.IsPrimeBigInt(p) {
		return nil, fmt.Errorf("Given p %s for prime number as modulo of the curve is not prime.\n", p)
	}

	modA := new(big.Int).Mod(a, p)
	modB := new(big.Int).Mod(b, p)
	modP := new(big.Int).Set(p)

	return &EllipticCurve{modA, modB, modP}, nil
}

// IsNonSingular True iff the elliptic curve isn't singular i.e. not crossing itself and not pointed,
// reducing risks of cryptographic vulnerability.
func (ec *EllipticCurve) IsNonSingular() bool {
	// 4 * a ** 3
	four := big.NewInt(4)
	aCubed := new(big.Int).Exp(ec.a, big.NewInt(3), nil)
	q := new(big.Int).Mul(four, aCubed)

	// 27 * b ** 2
	twentySeven := big.NewInt(27)
	bSquared := new(big.Int).Exp(ec.b, big.NewInt(2), nil)
	d := new(big.Int).Mul(twentySeven, bSquared)

	q.Add(q, d) // 4 * a**3  +  27 * b**2

	q.Mod(q, ec.p) // q = q % p

	return q.Cmp(big.NewInt(0)) != 0
}

// CalcSlope Calculates the slope existing between the segment linking both given points.
// Takes into consideration that P could be Q.
func (ec *EllipticCurve) CalcSlope(p *Point, q *Point) *big.Int {
	if p.Equals(q) {
		// (3 * (p.x^2)  + a)
		three := big.NewInt(3)                                 // 3
		pxSquared := new(big.Int).Exp(p.x, big.NewInt(2), nil) // p.x^2
		left := new(big.Int).Mul(three, pxSquared)             // 3 * p.x^2
		left.Add(left, ec.a)                                   // 3 * p.x^2 + a
		left.Mod(left, ec.p)                                   // .. mod p

		// modular_inverse of 2 * p.y % p
		twoPy := new(big.Int).Mul(big.NewInt(2), p.y) // 2 * p.y
		right := new(big.Int).ModInverse(twoPy, ec.p) // mod_inverse of 2 * p.y mod p

		//in case of no modular_inverse (rare)
		if right == nil {
			return nil
		}

		left.Mul(left, right) // (3 * (p.x)^2 + a) * ( mod_inv(2 * p.y) mod p)
		left.Mod(left, ec.p)  // ^^^^^^ mod p
		return left

	} else {
		// P != Q

		//left
		qySubPy := new(big.Int).Sub(q.y, p.y) // q.y - p.y

		//right
		// modular_inverse of q.x - p.x mod ec.p
		qxSubPx := new(big.Int).Sub(q.x, p.x)
		right := new(big.Int).ModInverse(qxSubPx, ec.p) // mod_inverse of (q.x - p.x) mod ec.p

		//if no modular inverse (rare)
		if right == nil {
			return nil
		}

		qySubPy.Mul(qySubPy, right) // (q.y - p.y) * mod_inverse of (q.x - p.x) mod ec.p
		qySubPy.Mod(qySubPy, ec.p)  // ^^^^^^ mod ec.p

		return qySubPy
	}
}

// ResolveX Resolves the x coordinate of the R point (to get the result of the sum of P & Q points on the curve).
// p & q points must be on the curve, else it will fail or give unexpected results (e.g. a point not on the curve).
// If one of p or q are NOT on the curve, will return -1 as a coordinate, and an error.
func (ec *EllipticCurve) ResolveX(slope *big.Int, p *Point, q *Point) (*big.Int, error) {
	//checks for points being on the curve
	if !ec.PointIsOnCurve(p) {
		return big.NewInt(-1), fmt.Errorf("P point (%s, %s) is not on curve.\n", p.x, p.y)
	}
	if !ec.PointIsOnCurve(q) {
		return big.NewInt(-1), fmt.Errorf("Q point (%s, %s) is not on curve.\n", q.x, q.y)
	}

	//(slope^2 - p.x - q.x) % ec.p
	RX := new(big.Int).Exp(slope, big.NewInt(2), nil) // slope ^ 2
	RX.Sub(RX, p.x)                                   // slope ^ 2 - p.x
	RX.Sub(RX, q.x)                                   // slope ^ 2 - p.x - q.x
	RX.Mod(RX, ec.p)                                  // ... mod p
	return RX, nil
}

// ResolveY Resolves the y coordinate of the R point (to get the result of the sum of P & Q points on the curve).
// p & q points must be on the curve, else it will fail or give unexpected results (e.g. a point not on the curve).
// If p is NOT on the curve, will return -1 as a coordinate, and an error.
func (ec *EllipticCurve) ResolveY(slope *big.Int, p *Point, RX *big.Int) (*big.Int, error) {
	//check for point being on the curve
	if !ec.PointIsOnCurve(p) {
		return big.NewInt(-1), fmt.Errorf("P point (%s, %s) is not on curve.\n", p.x, p.y)
	}

	// (slope*(p.x-xR)-p.y) % ec.p
	YR := new(big.Int).Sub(p.x, RX) // (p.x - xR)
	YR.Mul(YR, slope)               // slope * (p.x - xR)
	YR.Sub(YR, p.y)                 // slope * (p.x - xR) - p.y
	YR.Mod(YR, ec.p)                // .. mod p

	return YR, nil
}

// PointIsOnCurve True iff the given point is on the elliptic curve, using its coordinates.
// i.e. this means that the point must fulfill the equation : y² = x³+ ax + b for the current curve.
func (ec *EllipticCurve) PointIsOnCurve(p *Point) bool {
	//if the point is the neutral element (we consider it to be a nil point)
	if p == nil {
		return true
	}

	//base equation y² = x³+ ax + b
	// y² mod p
	pYSquaredModP := new(big.Int).Exp(p.y, big.NewInt(2), ec.p) // y² mod p

	//x³ + ax + b
	xCube := new(big.Int).Exp(p.x, big.NewInt(3), ec.p) // x³ mod p
	AX := new(big.Int).Mul(ec.a, p.x)                   // a * p.x
	AX.Mod(AX, ec.p)                                    // a * p.x mod p

	res := new(big.Int).Add(xCube, AX) // x³ + a * x
	res.Add(res, ec.b)                 // x³ + a * x + b
	res.Mod(res, ec.p)                 // .. mod p

	//true iff y² = x³ + ax + b
	return pYSquaredModP.Cmp(res) == 0
}

// SumPointsOnCurve Returns the sum of the two points on the curve.
// Will check if both points are on the curve at start, and will return a nil pointer with an error message.
// It will also check that the elliptic curve is non-singular, and that the resulting point is also on the curve.
// If not, will return a nil pointer, and an error message.
// WARNING: if both nil as a point is returned, and same for error, this means that
// p & q points are aligned, therefore the resulting point is the omega (infinity)
func (ec *EllipticCurve) SumPointsOnCurve(p *Point, q *Point) (*Point, error) {
	if !ec.IsNonSingular() {
		return nil, fmt.Errorf("Elliptic curve was singular, hence the sum of p(%s,%s) and q(%s,%s) won't be done.\n", p.x, p.y, q.x, q.y)
	}

	//checks if both points are on the curve
	if !ec.PointIsOnCurve(p) {
		return nil, fmt.Errorf("P point (%s, %s) is not on curve.\n", p.x, p.y)
	}
	if !ec.PointIsOnCurve(q) {
		return nil, fmt.Errorf("Q point (%s, %s) is not on curve.\n", q.x, q.y)
	}

	//when p or q is infinity, since P + omega = P (or Q + omega = Q)
	if p == nil {
		return q, nil // p = nil (i.e. omega), so p + q <==> omega + q
	}

	if q == nil {
		return p, nil //q = nil (i.e. omega), so p + q <==> p + omega
	}

	// p.y + q.y mod p
	PYplusQY := new(big.Int).Add(p.y, q.y)
	PYplusQY.Mod(PYplusQY, ec.p)

	//checks if p & q are not aligned (if they are, this means that the segment between them is omega itself)
	if p.x.Cmp(q.x) == 0 && PYplusQY.Cmp(big.NewInt(0)) == 0 {
		return nil, nil
	}

	slope := ec.CalcSlope(p, q)
	XR, err := ec.ResolveX(slope, p, q)
	if err != nil {
		return nil, err
	}

	YR, err := ec.ResolveY(slope, p, XR)
	if err != nil {
		return nil, err
	}

	R, err := NewPoint(XR, YR)
	if err != nil {
		return nil, err
	}

	return R, nil
}

// MultiplyPointByScalar Returns the point of the curve resulting
// of the multiplication of the given point of the curve by a scalar.
// Uses the double and add algorithm by using binary exponentiation.
// WARNING: if multiplying by zero, will return nil (as the omega) and nil (for error).
// WARNING: the parameter n will be modified.
func (ec *EllipticCurve) MultiplyPointByScalar(p *Point, n *big.Int) (*Point, error) {
	if !ec.PointIsOnCurve(p) {
		return nil, fmt.Errorf("Given point p(%s,%s) is not on the curve. Aborting multiplication.\n", p.x, p.y)
	}

	//if multiplying the point by zero, return omega (i.e. nil)
	if n.Cmp(big.NewInt(0)) == 0 {
		return nil, nil
	}

	//starting with the omega (neutral element)
	var res *Point = nil
	current := p.CopyPoint()

	zero := big.NewInt(0)
	two := big.NewInt(2)
	one := big.NewInt(1)

	//while n > 0
	for n.Cmp(zero) == 1 {

		modCond := new(big.Int).Mod(n, two) // n mod 2
		// if n mod 2 == 1
		if modCond.Cmp(one) == 0 {
			// res + current (first occurrence will be omega + current)
			tmp, err := ec.SumPointsOnCurve(res, current)
			if err != nil {
				fmt.Printf("Error during SumPointsOnCurve with res = %v and current = %v: %v", res, current, err)
				return nil, err
			}
			res = tmp
		}
		//current + current
		tmp, err := ec.SumPointsOnCurve(current, current)
		if err != nil {
			fmt.Printf("Error during SumPointsOnCurve with current = %v and current = %v: %v", current, current, err)
			return nil, err
		}
		current = tmp
		n.Div(n, two)
	}
	return res, nil
}

// CheckHasseTheorem True iff the elliptic curve passes the hasse theorem stating that the number of points of the curve
//
//	within its finite field is close to the number of elements of the finite field itself
//
// i.e. |N - (q + 1)| <= 2 * sqrt(p), see https://en.wikipedia.org/wiki/Hasse%27s_theorem_on_elliptic_curves
func (ec *EllipticCurve) CheckHasseTheorem() bool {
	//TODO: use schoof to get the number of points of the elliptic curve
	return ec.IsNonSingular() && true
}

func CreateEC() *EllipticCurve {
	a, ok := new(big.Int).SetString("0x0000000000000000000000000000000000000000000000000000000000000000", 0)
	if !ok {
		log.Fatalf("Given 'a' hexadecimal string representation was not correct or could not be converted to big.Int type.\n")
	}
	log.Printf("Successfully loaded 'a': %x \n", a)

	b, ok := new(big.Int).SetString("0x0000000000000000000000000000000000000000000000000000000000000007", 0)
	if !ok {
		log.Fatalf("Given 'b' hexadecimal string representation was not correct or could not be converted to big.Int type.\n")
	}
	log.Printf("Successfully loaded 'b': %x \n", b)

	p, ok := new(big.Int).SetString("0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F", 0)
	if !ok {
		log.Fatalf("Given 'p' hexadecimal string representation was not correct or could not be converted to big.Int type.\n")
	}
	log.Printf("Successfully loaded 'p': %x \n", p)

	curve, _ := NewEllipticCurve(a, b, p)
	return curve
}

// ProcessYFrom Computes Weiestrass equation to get Y
// if no square root, returns nil with false
func (ec *EllipticCurve) ProcessYFrom(x *big.Int) (*big.Int, bool) {
	x3 := new(big.Int).Exp(x, big.NewInt(3), ec.p) // x^3 mod p
	ax := new(big.Int).Mul(ec.a, x)                //a * x
	ax.Mod(ax, ec.GetP())                          // a * x mod p

	val := new(big.Int).Add(x3, ax)      // x^3 + a*x
	val.Add(val, ec.b)                   //x^3 + a*x + b
	val.Mod(val, ec.p)                   // .. mod p
	y := new(big.Int).ModSqrt(val, ec.p) // sqrt(val) mod p

	if y == nil {
		return nil, false
	}
	return y, true
}
