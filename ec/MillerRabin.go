package ec

import (
	"fmt"
	"math/big"
	"math/rand"
	"time"
)

func MillerRabin(candidate *big.Int) bool {
	// Simple shortcuts.
	one := big.NewInt(1)
	two := big.NewInt(2)

	// Read the candidate from the first argument, or default to 221 for test
	// purposes.

	modulo := new(big.Int)
	modulo.Sub(candidate, one)

	// Write the modulo (candidate -1) number in the form
	// 2^s * d.
	s := big.NewInt(0)
	remainder := new(big.Int)
	quotient := new(big.Int)
	quotient.Set(modulo)

	for remainder.Sign() == 0 {
		quotient.DivMod(quotient, two, remainder)
		s.Add(s, one)
	}
	// The last division failed, so we must decrement `s`.
	s.Sub(s, one)
	// quotient here contains the leftover which we could not divide by two,
	// and we have a 1 remaining from this last division.
	d := big.NewInt(1)
	d.Add(one, d.Mul(two, quotient))

	// Random number source for generating witnesses.
	r := rand.New(rand.NewSource(time.Now().UnixNano()))

	// Here 10 is the precision. Every increment to this value decreases the
	// chance of a false positive by 3/4.
	for k := 0; k < 10; k++ {

		// Every witness may prove that the candidate is composite, or assert
		// nothing.
		witness := new(big.Int)
		witness.Rand(r, modulo)

		exp := new(big.Int)
		exp.Set(d)
		generated := new(big.Int)
		generated.Exp(witness, exp, candidate)

		if generated.Cmp(modulo) == 0 || generated.Cmp(one) == 0 {
			continue
		}

		s64 := s.Int64()
		sInt := int(s64)

		for i := 1; i < sInt; i++ {
			generated.Exp(generated, two, candidate)

			if generated.Cmp(one) == 0 {
				fmt.Println("MillerRabin: Composite, not prime.")
				return false
			}

			if generated.Cmp(modulo) == 0 {
				break
			}
		}

		if generated.Cmp(modulo) != 0 {
			// We arrived here because the `i` loop ran its course naturally
			// without meeting the `x == modulo` break.
			fmt.Println("MillerRabin: Composite, not prime.")
			return false
		}
	}
	return true
}
