import sage.all

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer

import sage.arith.misc as arith

import webbrowser
from functools import reduce
from typing import List
from itertools import islice


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
