package utils

import (
	"math/big"
)

func IsPrime(n int64) bool {
	bn := big.NewInt(n)
	return bn.ProbablyPrime(20)
}
func IsPrimeBigInt(n *big.Int) bool {
	return n.ProbablyPrime(20)
}

func CreateSmallPrime(n int64) []*big.Int {
	var primes []*big.Int
	for i := int64(2); i < n; i++ {
		if IsPrime(i) {
			primes = append(primes, big.NewInt(i))
		}
	}
	return primes
}
