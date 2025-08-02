package main

import (
	"goschoof/ec"
	"log"
	"math/big"
)

func main() {

	//For testing purposes only

	log.SetPrefix("ec -- ")

	//secp256k1 curve
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

	curve, err := ec.NewEllipticCurve(a, b, p)
	if err != nil {
		log.Fatalf(err.Error())
	}

	log.Printf("Successfully created 'curve' EllipticCurve: (addresses) %x \n", curve)

	log.Printf("Checking that the neutral point (nil) is on the curve: %t\n", curve.PointIsOnCurve(nil))

	//secp256k1 generator point, see https://en.bitcoin.it/wiki/Secp256k1
	x, ok := new(big.Int).SetString("0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798", 0)
	y, ok := new(big.Int).SetString("0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8", 0)
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
}
