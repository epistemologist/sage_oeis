import sage.all

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer

import sage.arith.misc as arith

import webbrowser
from functools import reduce
from math import inf
from typing import List
from itertools import islice, count


class OEISSequence(SageObject):
    def __init__(self,
                 seq_number: Integer,
                 description: str,
                 offset: Integer = 1,
                 all_at_once: bool = True
                 ):
        self.seq_number = Integer(seq_number)
        self.description = description
        self.offset = offset
        self.url = f"http://oeis.org/A{self.seq_number}"
        # whether or not to generate all terms up to n or just term given
        self.all_at_once = all_at_once

    def __repr__(self):
        return self.description

    def __eq__(self, other):
        return isinstance(other, OEISSequence) and self.seq_number == other.seq_number

    def __ne__(self, other):
        return not self.__eq__(other)

    def _sage_src_(self):
        from sage.misc.sageinspect import sage_getsource
        return sage_getsource(self.__class__)

    # For evlauating a single sequence element
    def _eval(self, n: Integer) -> Integer:
        if self.all_at_once:
            return self._eval_up_to_n(n)[-1]
        else:
            raise NotImplementedError

    # For evaluating all sequence elements up to n (may help with speed)
    def _eval_up_to_n(self, n: Integer) -> List[Integer]:
        if not self.all_at_once:
            return [self._eval(i) for i in range(self.offset, self.offset+n)]
        else:
            raise NotImplementedError

    def __call__(self, n):
        return self._eval(n)

    def open_in_browser(self):
        webbrowser.open(self.url)


# Example sequence
class A000027(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            seq_number=27,
            description="The positive integers. Also called the natural numbers, the whole numbers or the counting numbers, but these terms are ambiguous.",
            all_at_once=False
        )

    def _eval(self, n) -> Integer:
        return n


class A000001(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1,
            description="Number of groups of order n.",
            all_at_once=False
        )

    def _eval(self, n) -> Integer:
        # Code originally from Jaap Spies (2007-02-04)
        from sage.libs.gap.libgap import libgap
        return 0 if n == 0 else Integer(libgap.NrSmallGroups(n))


class A000002(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            seq_number=2,
            description="Kolakoski sequence: a(n) is length of n-th run; a(1) = 1; sequence consists just of 1's and 2's.",
            all_at_once=True
        )

    def _kolakoski_iterator(self):
        # From https://11011110.github.io/blog/2016/10/14/kolakoski-sequence-via.html
        x = y = -1
        while True:
            yield [2, 1][x & 1]
            f = y & ~ (y+1)
            x ^= f
            y = (y+1) | (f & (x >> 1))

    def _eval_up_to_n(self, n: Integer) -> List[Integer]:
        return [Integer(i) for i in islice(self._kolakoski_iterator(), n)]


class A000004(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            seq_number=4,
            description="The zero sequence.",
            all_at_once=False
        )

    def _eval(self, n) -> Integer:
        return Integer(0)


class A000005(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            seq_number=5,
            description="d(n) (also called tau(n) or sigma_0(n)), the number of divisors of n.",
            all_at_once=False
        )

    def _eval(self, n) -> Integer:
        return arith.number_of_divisors(n)


class A000007(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=7,
            description="The characteristic function of {0}: a(n) = 0^n.",
            all_at_once=False
        )

    def _eval(self, n) -> Integer:
        return Integer(1) if n == 0 else Integer(0)


class A000009(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=9,
            description="Expansion of Product_{m >= 1} (1 + x^m); number of partitions of n into distinct parts; number of partitions of n into odd parts.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Slow way to compute terms - there is faster recurrence on the OEIS page
        from sage.all import PolynomialRing, ZZ
        R = PolynomialRing(ZZ, names=('x',))
        (x,) = R.gens()
        product = reduce(lambda x, y: x*y, [1 + x**i for i in range(1, n+1)])
        return product.coefficients()[::-1][:n]


class A000010(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=10,
            description="Euler totient function phi(n): count numbers <= n and prime to n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return arith.euler_phi(n)


class A000012(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=12,
            description="The simplest sequence of positive numbers: the all 1's sequence.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(1)

# TODO: A000014


class A000019(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=19,
            description="Number of primitive permutation groups of degree n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.libs.gap.libgap import libgap
        return 1 if n == 1 else Integer(libgap.NrPrimitiveGroups(n))


class A000029(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=29,
            description="Number of necklaces with n beads of 2 colors, allowing turning over (these are also called bracelets).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        if n == 0:
            return 1
        else:
            return sum([
                (arith.euler_phi(d) * pow(2, n/d)) / (2*n)
                for d in arith.divisors(n)
            ]) + (pow(2, (n-1)//2) if n % 2 == 1 else (pow(2, n//2-1) + pow(2, n//2 - 2)))


class A000031(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=31,
            description="Number of n-bead necklaces with 2 colors when turning over is not allowed; also number of output sequences from a simple n-stage cycling shift register; also number of binary irreducible polynomials whose degree divides n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        if n == 0:
            return 1
        else:
            return sum([arith.euler_phi(d) * pow(2, n/d) for d in arith.divisors(n)]) // n


class A000032(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=32,
            description="Lucas numbers beginning at 2: L(n) = L(n-1) + L(n-2), L(0) = 2, L(1) = 1.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List[Integer]:
        seq = [2, 1]
        while len(seq) < n:
            seq.append(seq[-1] + seq[-2])
        return seq[:n]


class A000035(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=35,
            description="Period 2: repeat [0, 1]; a(n) = n mod 2; parity of n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return n % 2


class A000040(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=40,
            description="The prime numbers.",
            all_at_once=True,
        )

    def _eval_up_to_n(self, n: Integer) -> List[Integer]:
        return arith.primes_first_n(n)


class A000041(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            seq_number=41,
            offset=0,
            description="a(n) is the number of partitions of n (the partition numbers).",
            all_at_once=False,
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.combinat.partition import Partitions
        return Partitions(n).cardinality()


class A000043(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            seq_number=43,
            offset=1,
            description="Mersenne exponents: primes p such that 2^p - 1 is prime. Then 2^p - 1 is called a Mersenne prime.",
            all_at_once=True
        )
        self.cached = cached
        self.SEQ_TERMS = [2, 3, 5, 7, 13, 17, 19, 31, 61, 89, 107, 127, 521, 607, 1279, 2203, 2281, 3217, 4253, 4423, 9689, 9941, 11213, 19937, 21701, 23209, 44497, 86243, 110503,
                          132049, 216091, 756839, 859433, 1257787, 1398269, 2976221, 3021377, 6972593, 13466917, 20996011, 24036583, 25964951, 30402457, 32582657, 37156667, 42643801, 43112609, 57885161]

    def _eval_up_to_n(self, n: Integer) -> List[Integer]:
        if self.cached:
            return [Integer(i) for i in self.SEQ_TERMS[:n]]
        else:
            seq = []
            for p in arith.primes(1, inf):
                if arith.is_prime(pow(2, p) - 1):
                    seq.append(p)
                if len(seq) == n:
                    return seq


class A000045(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=45,
            description="Fibonacci numbers: F(n) = F(n-1) + F(n-2) with F(0) = 0 and F(1) = 1.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.combinat.combinat import fibonacci
        return fibonacci(n)


class A000048(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=48,
            description="Number of n-bead necklaces with beads of 2 colors and primitive period n, when turning over is not allowed but the two colors can be interchanged.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        if n == 0:
            return 1
        else:
            return sum([arith.moebius(d)*pow(2, n//d) for d in arith.divisors(n) if d % 2 == 1])//(2*n)


# Some sequences dealing with the enumeration of trees ()
class A000081(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=81,
            description="Number of unlabeled rooted trees with n nodes (or connected functions with a fixed point).",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List[Integer]:
        # Adapted from "The Art of Computer Programming" Vol 1 2.3.4.4 Exercise 2
        a = [0,1]
        s = lambda n,k: sum([a[n+1-j*k] for j in range(1, n//k+1)])
        for n_ in range(1, n):
            a.append(sum([i * a[i] * s(n_, i) for i in range(1, n_+1)]) // n_)
        return [Integer(i) for i in a][:n]

class A000055(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=55,
            description="Number of trees with n unlabeled nodes.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List[Integer]:
        # Use the generating function relation 
        # A(x) = 1 + T(x) - T(x)^2 / 2 + T(x^2) / 2, where T(x) is the g.f. of A000081 
        terms = A000081()._eval_up_to_n(50)
        from sage.all import PolynomialRing, ZZ
        R = PolynomialRing(ZZ, names=('x',))
        (x,) = R.gens()
        T = sum([a*x**i for i, a in enumerate(terms) ] )
        A = 1 + T - (T*T)/2 + T.subs({x: x*x})/2
        return [Integer(i) for i in A.coefficients()[:n]]
    
class A000058(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=58,
            description="Sylvester's sequence: a(n+1) = a(n)^2 - a(n) + 1, with a(0) = 2.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        seq = [2]
        while len(seq) < n:
            seq.append(seq[-1]**2 - seq[-1] + 1)
        return [Integer(i) for i in seq[:n]]

class A000069(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=69,
            description="Odious numbers: numbers with an odd number of 1's in their binary expansion.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count(1):
            if bin(i).count('1') % 2 == 1:
                seq.append(i)
            if len(seq) >= n: break
        return [Integer(i) for i in seq][:n]
    
class A000079(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            seq_number=79,
            offset=0,
            description="Powers of 2: a(n) = 2^n.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(1 << n)
