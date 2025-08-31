package polynom

import (
	"fmt"
	"math/big"
)

type Polynom struct {
	Coefficients []*big.Int
	P            *big.Int
}

func NewPolynom(coeffs []*big.Int, p *big.Int) *Polynom {
	c := make([]*big.Int, len(coeffs))
	// put all the coefficients mod p for finite field integration
	for i, coeff := range coeffs {
		c[i] = new(big.Int).Mod(coeff, p)
	}
	return &Polynom{
		Coefficients: c,
		P:            new(big.Int).Set(p),
	}
}

// Degree returns the degree of the polynom (last non-null coefficient)
func (poly *Polynom) Degree() int {
	for i := len(poly.Coefficients) - 1; i >= 0; i-- {
		if poly.Coefficients[i].Sign() != 0 {
			return i
		}
	}
	// polynom nul = degree 0 by convention
	return 0
}

func (poly *Polynom) String() string {
	deg := poly.Degree()
	s := ""
	for i := deg; i >= 0; i-- {
		coeff := poly.Coefficients[i]
		if coeff.Sign() == 0 {
			continue
		}
		if s != "" {
			s += " + "
		}
		s += fmt.Sprintf("%sx^%d", coeff.String(), i)
	}
	if s == "" {
		s = "0"
	}
	return s
}

// Add adds two polynomials mod p & return the resulting polynom.
func (poly *Polynom) Add(other *Polynom) *Polynom {
	maxLen := len(poly.Coefficients)
	if len(other.Coefficients) > maxLen {
		maxLen = len(other.Coefficients)
	}

	resultCoeffs := make([]*big.Int, maxLen)
	zero := big.NewInt(0)

	// add all the corresponding coefficients
	for i := 0; i < maxLen; i++ {
		a := zero
		b := zero
		if i < len(poly.Coefficients) {
			a = poly.Coefficients[i]
		}
		if i < len(other.Coefficients) {
			b = other.Coefficients[i]
		}

		sum := new(big.Int).Add(a, b)
		sum.Mod(sum, poly.P)
		resultCoeffs[i] = sum
	}
	return &Polynom{
		Coefficients: resultCoeffs,
		P:            new(big.Int).Set(poly.P),
	}
}

func (poly *Polynom) Sub(other *Polynom) *Polynom {
	maxLen := len(poly.Coefficients)
	if len(other.Coefficients) > maxLen {
		maxLen = len(other.Coefficients)
	}

	resultCoeffs := make([]*big.Int, maxLen)
	zero := big.NewInt(0)

	// subtract all the corresponding coefficients
	for i := 0; i < maxLen; i++ {
		a := zero
		b := zero
		if i < len(poly.Coefficients) {
			a = poly.Coefficients[i]
		}
		if i < len(other.Coefficients) {
			b = other.Coefficients[i]
		}

		diff := new(big.Int).Sub(a, b)
		diff.Mod(diff, poly.P)
		resultCoeffs[i] = diff
	}
	return &Polynom{
		Coefficients: resultCoeffs,
		P:            new(big.Int).Set(poly.P),
	}
}

func (poly *Polynom) Mul(other *Polynom) *Polynom {
	degA := len(poly.Coefficients)
	degB := len(other.Coefficients)
	resultDegree := degA + degB - 1

	resultCoeffs := make([]*big.Int, resultDegree)
	for i := range resultCoeffs {
		resultCoeffs[i] = big.NewInt(0)
	}

	for i := 0; i < degA; i++ {
		for j := 0; j < degB; j++ {
			prod := new(big.Int).Mul(poly.Coefficients[i], other.Coefficients[j])
			prod.Mod(prod, poly.P)
			resultCoeffs[i+j].Add(resultCoeffs[i+j], prod)
			resultCoeffs[i+j].Mod(resultCoeffs[i+j], poly.P)
		}
	}

	return &Polynom{
		Coefficients: resultCoeffs,
		P:            new(big.Int).Set(poly.P),
	}
}

func (poly *Polynom) Copy() *Polynom {
	if poly == nil {
		return nil
	}
	coeffs := make([]*big.Int, len(poly.Coefficients))
	for i, c := range poly.Coefficients {
		if c == nil {
			coeffs[i] = nil
		} else {
			coeffs[i] = new(big.Int).Set(c)
		}
	}
	pCopy := new(big.Int)
	if poly.P != nil {
		pCopy.Set(poly.P)
	}
	return &Polynom{
		Coefficients: coeffs,
		P:            pCopy,
	}
}

func (poly *Polynom) trimTrailingZeros() {
	i := len(poly.Coefficients) - 1
	for i > 0 && poly.Coefficients[i].Sign() == 0 {
		i--
	}
	poly.Coefficients = poly.Coefficients[:i+1]
}

func (poly *Polynom) IsZero() bool {
	for _, c := range poly.Coefficients {
		if c.Sign() != 0 {
			return false
		}
	}
	return true
}

func (poly *Polynom) NormalizeMonic() (*Polynom, *big.Int) {
	hn := poly.Copy()
	hn.trimTrailingZeros()
	deg := hn.Degree()
	lead := hn.Coefficients[deg]
	inv := new(big.Int).ModInverse(lead, hn.P)
	if inv == nil {
		return nil, nil
	}
	for i := 0; i < len(hn.Coefficients); i++ {
		if hn.Coefficients[i] == nil {
			hn.Coefficients[i] = big.NewInt(0)
		}
		hn.Coefficients[i].Mul(hn.Coefficients[i], inv)
		hn.Coefficients[i].Mod(hn.Coefficients[i], hn.P)
	}
	return hn, inv
}

func (poly *Polynom) DivMod(h *Polynom) (*Polynom, *Polynom) {
	if h == nil || h.IsZero() {
		return NewPolynom([]*big.Int{big.NewInt(0)}, poly.P), poly.Copy()
	}
	R := poly.Copy()
	R.trimTrailingZeros()
	Q := NewPolynom([]*big.Int{big.NewInt(0)}, poly.P)

	hMonic, _ := h.NormalizeMonic()
	hMonic.trimTrailingZeros()
	degH := hMonic.Degree()

	// if deg(F) < deg(h), q 0, r F
	if R.Degree() < degH {
		return Q, R
	}

	ensureLen := func(a *[]*big.Int, n int) {
		if len(*a) < n {
			old := *a
			extra := make([]*big.Int, n-len(old))
			for i := range extra {
				extra[i] = big.NewInt(0)
			}
			*a = append(old, extra...)
		}
	}

	p := new(big.Int).Set(poly.P)

	for R.Degree() >= degH && !R.IsZero() {
		degR := R.Degree()
		k := degR - degH
		c := new(big.Int).Set(R.Coefficients[degR])

		// Q[k] += c
		ensureLen(&Q.Coefficients, k+1)
		Q.Coefficients[k].Add(Q.Coefficients[k], c)
		Q.Coefficients[k].Mod(Q.Coefficients[k], p)

		// R -= c * x^k * hMonic
		for i := 0; i <= degH; i++ {
			idx := i + k
			if idx >= len(R.Coefficients) {
				extra := make([]*big.Int, idx-len(R.Coefficients)+1)
				for j := range extra {
					extra[j] = big.NewInt(0)
				}
				R.Coefficients = append(R.Coefficients, extra...)
			}
			if R.Coefficients[idx] == nil {
				R.Coefficients[idx] = big.NewInt(0)
			}
			tmp := new(big.Int).Mul(c, hMonic.Coefficients[i])
			tmp.Mod(tmp, p)
			R.Coefficients[idx].Sub(R.Coefficients[idx], tmp)
			R.Coefficients[idx].Mod(R.Coefficients[idx], p)
		}
		R.trimTrailingZeros()
	}

	Q.trimTrailingZeros()
	R.trimTrailingZeros()
	return Q, R
}

func (poly *Polynom) PowMod(n *big.Int, h *Polynom) *Polynom {
	one := NewPolynom([]*big.Int{big.NewInt(1)}, poly.P)
	if n.Sign() == 0 {
		return one
	}
	base := poly.Copy()
	_, base = base.DivMod(h)
	res := one
	e := new(big.Int).Set(n)
	zero := big.NewInt(0)
	two := big.NewInt(2)

	for e.Cmp(zero) > 0 {
		if new(big.Int).And(e, big.NewInt(1)).Cmp(zero) != 0 {
			res = res.Mul(base)
			_, res = res.DivMod(h)
		}
		base = base.Mul(base)
		_, base = base.DivMod(h)
		e.Div(e, two)
	}
	return res
}

func PolyInvMod(f, h *Polynom) (*Polynom, bool) {
	// Extended GCD: uf + vh = g
	g, u, _ := PolyExtGCD(f, h)
	if g.Degree() != 0 || !g.LeadingCoeffIsOne() {
		c := g.Coeff(0)
		if c.Sign() == 0 {
			return nil, false
		}
		cInv := new(big.Int).ModInverse(c, f.P)
		if cInv == nil {
			return nil, false
		}
		u = u.Scale(cInv)
		u = u.ModCoeffs() // reduce coeffs mod p
	}
	_, u = u.DivMod(h)
	return u, true
}

func GCDPolynom(a, b *Polynom) *Polynom {
	A := a.Copy()
	B := b.Copy()
	for !B.IsZero() {
		_, R := A.DivMod(B)
		A, B = B, R
	}

	lc := A.LeadingCoeff()
	if lc == nil || lc.Sign() == 0 {
		return A
	}
	lcInv := new(big.Int).ModInverse(lc, A.P)
	if lcInv != nil {
		A = A.Scale(lcInv)
		A = A.ModCoeffs()
	}
	return A
}

func PolyExtGCD(a, b *Polynom) (*Polynom, *Polynom, *Polynom) {
	// x2=1, x1=0; y2=0, y1=1
	one := NewPolynom([]*big.Int{big.NewInt(1)}, a.P)
	zero := NewPolynom([]*big.Int{big.NewInt(0)}, a.P)

	r2, r1 := a.Copy(), b.Copy()
	x2, x1 := one.Copy(), zero.Copy()
	y2, y1 := zero.Copy(), one.Copy()

	for !r1.IsZero() {
		q, r := r2.DivMod(r1)
		r2, r1 = r1, r
		// x = x2 - q*x1
		x := x2.Sub(q.Mul(x1))
		x = x.ModCoeffs()
		// y = y2 - q*y1
		y := y2.Sub(q.Mul(y1))
		y = y.ModCoeffs()
		x2, x1 = x1, x
		y2, y1 = y1, y
	}
	// normalise g = r2
	lc := r2.LeadingCoeff()
	if lc != nil && lc.Sign() != 0 {
		lcInv := new(big.Int).ModInverse(lc, a.P)
		if lcInv != nil {
			r2 = r2.Scale(lcInv).ModCoeffs()
			x2 = x2.Scale(lcInv).ModCoeffs()
			y2 = y2.Scale(lcInv).ModCoeffs()
		}
	}
	return r2, x2, y2
}

func DivExact(a, b *Polynom) *Polynom {
	q, r := a.DivMod(b)
	if !r.IsZero() {
		panic("DivExact: remainder not zero")
	}
	return q
}

func (poly *Polynom) LeadingCoeff() *big.Int {
	if poly.Degree() < 0 {
		return big.NewInt(0)
	}
	return new(big.Int).Set(poly.Coeff(poly.Degree()))
}

// test if dominant coeff is 1 mod p
func (poly *Polynom) LeadingCoeffIsOne() bool {
	lc := poly.LeadingCoeff()
	if lc == nil {
		return false
	}
	m := new(big.Int).Mod(lc, poly.P)
	return m.Cmp(big.NewInt(1)) == 0
}

// reduces all coefficients mod p
func (poly *Polynom) ModCoeffs() *Polynom {
	for i := 0; i <= poly.Degree(); i++ {
		ci := new(big.Int).Mod(poly.Coeff(i), poly.P)
		if ci.Sign() < 0 {
			ci.Add(ci, poly.P)
		}
		poly.SetCoeff(i, ci)
	}
	poly.trimTrailingZeros()
	return poly
}

func (poly *Polynom) Coeff(i int) *big.Int {
	if i < 0 || i >= len(poly.Coefficients) {
		return big.NewInt(0)
	}
	return new(big.Int).Set(poly.Coefficients[i])
}

func (poly *Polynom) SetCoeff(i int, v *big.Int) {
	if i < 0 {
		return
	}
	if i >= len(poly.Coefficients) {
		extra := make([]*big.Int, i+1-len(poly.Coefficients))
		for k := range extra {
			extra[k] = big.NewInt(0)
		}
		poly.Coefficients = append(poly.Coefficients, extra...)
	}
	// v mod p
	vv := new(big.Int).Mod(v, poly.P)
	if vv.Sign() < 0 {
		vv.Add(vv, poly.P)
	}
	poly.Coefficients[i] = vv
}

func (poly *Polynom) Scale(k *big.Int) *Polynom {
	if poly == nil {
		return poly
	}
	km := new(big.Int).Mod(k, poly.P)
	if km.Sign() < 0 {
		km.Add(km, poly.P)
	}
	for i := 0; i < len(poly.Coefficients); i++ {
		if poly.Coefficients[i] == nil {
			poly.Coefficients[i] = big.NewInt(0)
			continue
		}
		poly.Coefficients[i].Mul(poly.Coefficients[i], km)
		poly.Coefficients[i].Mod(poly.Coefficients[i], poly.P)
		if poly.Coefficients[i].Sign() < 0 {
			poly.Coefficients[i].Add(poly.Coefficients[i], poly.P)
		}
	}
	poly.trimTrailingZeros()
	return poly
}
