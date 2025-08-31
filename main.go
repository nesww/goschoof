package main

import (
	"goschoof/ec"
	"goschoof/polynom"
	"goschoof/schoof"
	"goschoof/utils"
	"log"
	"math/big"
)

func main() {

	//For testing purposes only

	log.SetPrefix("ec -- ")

	//secp256k1 curve
	curve := ec.CreateEC()

	log.Printf("Successfully created 'curve' EllipticCurve: (addresses) %x \n", curve)

	log.Printf("Checking that the neutral point (nil) is on the curve: %t\n", curve.PointIsOnCurve(nil))

	//secp256k1 generator point, see https://en.bitcoin.it/wiki/Secp256k1
	x, _ := new(big.Int).SetString("0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 0)
	y, _ := new(big.Int).SetString("0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 0)
	gen, err := ec.NewPoint(x, y)

	log.Printf("Checking that the point 'gen'\n(%s, \n%s)\n is on the curve: %t\n", x, y, curve.PointIsOnCurve(gen))

	zero := new(big.Int).Set(big.NewInt(0))
	z, err := ec.NewPoint(zero, zero)

	log.Printf("Checking that the point (0,0) is on the curve: %t\n", curve.PointIsOnCurve(z))

	res, err := curve.MultiplyPointByScalar(gen, big.NewInt(2))
	if err != nil {
		log.Fatalf(err.Error())
	}
	log.Printf("Result of the addition of multiplying our 'gen' by 2: \n(%s,\n%s)\n", res.GetX(), res.GetY())

	////////////////////////////////////
	// Frobenius					  //
	////////////////////////////////////

	log.Printf("-------------------------------------------\n")
	log.Printf("Starting to generate the prime numbers for the curve and testing Frobenius\n")

	primes := utils.CreateSmallPrime(20)
	log.Printf("Prime numbers: %v", primes)

	test, _ := ec.NewPoint(big.NewInt(21), big.NewInt(82))
	ta, tb := test.Frobenius(big.NewInt(19))
	log.Printf("Frobenius results for test(%v,%v): %v, %v", test.GetX(), test.GetY(), ta, tb)

	a := big.NewInt(21)
	b := big.NewInt(7689)
	p := big.NewInt(83)

	curve2, err := ec.NewEllipticCurve(a, b, p)
	if err != nil {
		log.Fatalf(err.Error())
	}

	points := schoof.CountTorsion2PointsFromPoly(curve2)

	log.Printf("Points of order 2 on the curve: %d found\n", points)

	points3 := schoof.ResolvePolynomialDivisionL3(curve2)

	log.Printf("Points of order 3 on the curve: %d found\n", len(points3))
	for i, pt := range points3 {
		log.Printf("Point %d : x = %s, y = %s\n", i+1, pt.GetX(), pt.GetY())
	}

	//using a cache to prevent stack overflow and better performance
	cache := schoof.NewPSICache(curve)
	l := big.NewInt(3)
	poly := schoof.PSI_l(curve, l, cache)
	log.Printf("Polynom found for l=%d %v", l, poly)

	//testing mod div for polynomials
	ptest := big.NewInt(101)
	poly = polynom.NewPolynom([]*big.Int{}, ptest)
	h := schoof.PSI_l(curve, l, cache)

	//x^3 + 2x + 5
	F := polynom.NewPolynom([]*big.Int{
		big.NewInt(5), // x^0
		big.NewInt(2), // x^1
		big.NewInt(0), // x^2
		big.NewInt(1), // x^3
	}, ptest)

	Q, R := F.DivMod(h)
	log.Printf("Degree(R) < Degree(h)?: %v", R.Degree() < h.Degree())

	log.Printf("R = %v", R)

	F2 := Q.Mul(h).Add(R) // Q*h + R (quotient * h + remainder must give the F polynomial)
	log.Printf("F2 == F ? \n %v == %v", F2, F)

	p = big.NewInt(101)
	h = polynom.NewPolynom([]*big.Int{big.NewInt(1), big.NewInt(0), big.NewInt(1)}, p) // x^2+1

	Xplus3 := polynom.NewPolynom([]*big.Int{big.NewInt(3), big.NewInt(1)}, p)     // x+3
	FourxPlus7 := polynom.NewPolynom([]*big.Int{big.NewInt(7), big.NewInt(4)}, p) // 4x+7

	F = h.Mul(Xplus3).Add(FourxPlus7)
	Q, R = F.DivMod(h)
	log.Printf("R expected ~ 4x+7, got R= %s\n\n", R.String())

	//////////////////
	// SCHOOF 		//
	//////////////////

	N := schoof.Schoof(curve2)
	log.Printf("N found for secp256k1: %v", N)
}
