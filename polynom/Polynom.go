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
	//put all the coefficients mod p for finite field integration
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
	//polynom nul = degree 0 by convention
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

// Add adds two polynoms mod p & return the resulting polynom.
func (poly *Polynom) Add(other *Polynom) *Polynom {
	maxLen := len(poly.Coefficients)
	if len(other.Coefficients) > maxLen {
		maxLen = len(other.Coefficients)
	}

	resultCoeffs := make([]*big.Int, maxLen)
	zero := big.NewInt(0)

	//add all the corresponding coefficients
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

	//subtract all the corresponding coefficients
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
	degA := poly.Degree()
	degB := other.Degree()
	resultDegree := degA + degB - 1

	resultCoeffs := make([]*big.Int, resultDegree)
	for i := range resultCoeffs {
		resultCoeffs[i] = big.NewInt(0)
	}

	for i := 0; i < degA; i++ {
		for j := 0; j < degB; i++ {
			prod := new(big.Int).Mul(poly.Coefficients[i], other.Coefficients[j])
			prod.Mod(prod, poly.P)
			resultCoeffs[i+j].Add(resultCoeffs[i+j], prod)
			resultCoeffs[i+j].Mod(resultCoeffs[i+j], poly.P)
		}
	}

	return &Polynom{
		Coefficients: resultCoeffs,
		P:            new(big.Int).Set(poly.P)}
}
