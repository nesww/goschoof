package ec

import (
	"fmt"
	"math/big"
)

type Point struct {
	x *big.Int
	y *big.Int
}

func NewPoint(x, y *big.Int) (*Point, error) {
	if x == nil || y == nil {
		return nil, fmt.Errorf("Coordinates cannot be nil.\n")
	}
	return &Point{x, y}, nil
}

func (p *Point) Equals(q *Point) bool {
	if p == nil && q == nil {
		return true
	}

	if p == nil || q == nil {
		return false
	}

	return p.x == q.x && p.y == q.y
}

// CopyPoint Creates a deep copy of the given point.
func (p *Point) CopyPoint() *Point {
	if p == nil {
		return nil
	}

	copy := new(Point)

	if p.x != nil {
		copy.x = new(big.Int).Set(p.x)
	}

	if p.y != nil {
		copy.y = new(big.Int).Set(p.y)
	}

	return copy
}

func (p *Point) Frobenius(mod *big.Int) (*big.Int, *big.Int) {
	x := new(big.Int).Exp(p.x, mod, mod)
	y := new(big.Int).Exp(p.y, mod, mod)
	return x, y
}

func (p *Point) GetX() *big.Int {
	return p.x
}

func (p *Point) GetY() *big.Int {
	return p.y
}
