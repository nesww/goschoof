package ecdh

import ( 
	"fmt"
)

type Operation func(a, b *Point) *Point

type Group struct{
	Elements []*Point 
	Operation Operation
	Neutral *Point
}

// Instantiate a new group instance.
// If fails, Group instance will be empty, and a error will be given.
func NewGroup(elements []*Point, op Operation, neutral *Point) (Group, error) {
	valid, reason := validate(&Group{elements, op, neutral});
	if !valid || len(reason) > 0 {
		return Group{}, fmt.Errorf("Group could not be created: invalid params: %s", reason);
	}	
	return Group{elements, op, neutral}, nil
}

//TODO: put in a dedicated module/utils package
// Checks wether a given slice of Points contains the given Point
// (considering that one point will be in the slice, iff 
// they have the same coordinates (x and y as big.Int).
func ContainsPointSlice(slice []*Point, val *Point) bool {
	for _, item := range slice {
		if item.Equals(val) {
			return true
		}
	}
	return false;
}

// Validates that the given group fulfills all the following conditions:
// - Has closure: the operation on elements of the group returns another element of the group
// - Has associativity: the order of the operation between 3 or more elements does not change the result
// - Has a neutral element: One element X, considered as neutral will act as following: a is an element, a op X = a
// - Each element has one inverse within the same group: for a given b element, -b exists such as b + (- b) = X (neutral element)
// If it failed to fulfill one or many conditions, each of them will be appended in a string that is returned with the boolean.
func validate(g *Group) (bool, string) {
	reason := "";
	closure := has_closure(g, &reason);
	associativity := has_associativity(g, &reason)
	neutral := has_neutral(g, &reason)
	inverses := has_inverses(g, &reason)
	return len(reason) == 0 && (closure && associativity && neutral && inverses), reason;
}

// Checks wether the given group fulfills the closure condition.
// This means that the operation of the group must give a result
// that is also an element of the group, else, the group isn't closed,
// and thus, does not fulfill this condition.
//
// i.e, for each elements a and b of the group, the element a op b is 
// an element of the group
func has_closure(g *Group, reason *string) bool {
	for _, a := range g.Elements {
		for _, b := range g.Elements {
			res := g.Operation(a, b);
			if !ContainsPointSlice(g.Elements, res) {
				*reason = fmt.Sprintf("Closure: Operation result of %s and %s is not in elements list of the group, result found: %s. The group does not fulfill closure condition.\n", a, b, res)
				return false;
			}	
		}
	}
	return true
}

// Checks if the given group fulfills the associativity conditon.
// This means that the order of the operation between 3 or more 
// elements does not change the result.
func has_associativity(g *Group, reason *string) bool {
	for _, a := range g.Elements {
		for _, b := range g.Elements {
			for _, c := range g.Elements {
				left_res := g.Operation(g.Operation(a,b),c)
				right_res := g.Operation(a, g.Operation(b,c))
				if !left_res.Equals(right_res) {
					*reason = fmt.Sprintf(*reason + "Associativity: the result of '(a op b) op c' was not equal to ' a op (b op c)': (%s op %s) op %s = %s != %s op (%s op %s) = %s\n", a, b, c, left_res, a, b, c, right_res)
					return false;
				}
			}
		}
	}
	return true;
}

// Checks if the given group fulfills the existence of a neutral element condition.
// This means that with one element X, considered as neutral, 
// it will act as following: a is an element, a op X = a 
func has_neutral(g *Group, reason *string) bool {
	omega := g.Neutral;
	for _, a := range g.Elements {
		aom := g.Operation(a, omega)
		oma := g.Operation(omega, a)
		if !aom.Equals(a) || !oma.Equals(a) {
			*reason = fmt.Sprintf(*reason + "Neutral: the a op omega or omega op  a result was not what it was expected: %s op %s != %s | %s op %s = %s, expected as a result for both: %s", a, omega, aom, omega, a, oma, a);
			return false
		}
	}
	return true;	
}

// Checks if the given group fulfills the existence of an inverse for each element in the group condition.
// This means that for a given b element, -b exists such as b + (- b) = X (neutral element) 
func has_inverses(g *Group, reason *string) bool {
	omega := g.Neutral
	for _, a := range g.Elements {
		
		has_inverse_iter := false
		for _, b := range g.Elements {
			has_inverse_iter = has_inverse_iter || (g.Operation(a,b).Equals(omega) && g.Operation(b,a).Equals(omega))
			if has_inverse_iter {break};
		}
		
		if !has_inverse_iter {
			*reason = fmt.Sprintf(*reason + "Inverse: the element %s has no inverse in the group.", a)
			return false;
		}
	}
	return true;
}
