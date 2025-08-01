package ecdh

import (
	"math/big"
)

type Point struct{
	x *big.Int
	y *big.Int
}

func NewPoint(x, y *big.Int) (Point,error) {
	return Point{x,y}, nil
}

func (p *Point) Equals(q *Point) bool {
	if p == nil && q == nil {
		return true;
	}

	if p == nil || q == nil {
		return false;
	}

	return p.x == q.x && p.y == q.y;
}

