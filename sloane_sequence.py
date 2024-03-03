import sage.all

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer

import sage.arith.misc as arith

import webbrowser
from functools import cache
from math import inf
from typing import List
from itertools import islice, count

"""
Sequences to implement:
 - A000014: Number of series-reduced trees with n nodes.
 - A000105: Number of free polyominoes (or square animals) with n cells.
 - A000109: Number of simplicial polyhedra with n vertices; simple planar graphs with n vertices and 3n-6 edges; maximal simple planar graphs with n vertices; planar triangulations with n vertices; triangulations of the sphere with n vertices; 3-connected cubic planar graphs on 2n-4 vertices.
 - A000609: Number of threshold functions of n or fewer variables.
 - A000612: Number of P-equivalence classes of switching functions of n or fewer variables, divided by 2.
 - A000798: Number of different quasi-orders (or topologies, or transitive digraphs) with n labeled elements.
 - A001034: Orders of noncyclic simple groups (without repetition).
"""


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
        return arith.prod([1 + x**i for i in range(1, n+1)]).coefficients()[::-1][:n]


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
        a = [0, 1]
        def s(n, k): return sum([a[n+1-j*k] for j in range(1, n//k+1)])
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
        T = sum([a*x**i for i, a in enumerate(terms)])
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
            if len(seq) >= n:
                break
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


class A000085(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            seq_number=85,
            offset=0,
            description="Number of self-inverse permutations on n letters, also known as involutions; number of standard Young tableaux with n cells.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        # Original code from Jaap Spies (2007-02-03)
        return sum(arith.factorial(n) // (arith.factorial(n-2*k) * (2**k) * arith.factorial(k)) for k in range(n//2+1))


class A000088(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            seq_number=88,
            offset=0,
            description="Number of graphs on n unlabeled nodes.",
            all_at_once=False
        )
        self.cached = cached
        self.SEQ_TERMS = [1, 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168,
                          1018997864, 165091172592, 50502031367952,
                          29054155657235488, 31426485969804308768,
                          64001015704527557894928,
                          245935864153532932683719776,
                          1787577725145611700547878190848,
                          24637809253125004524383007491432768]

    def _eval(self, n: Integer) -> Integer:
        if self.cached:
            return Integer(self.SEQ_TERMS[n])
        else:
            # Use nauty-geng to count number of graphs on n nodes
            # On my machine, this is impractical for any n > 10

            # Code modified from https://github.com/sagemath/sage/blob/develop/src/sage/graphs/graph_generators.py
            from sage.features.nauty import NautyExecutable
            import subprocess
            import re

            geng_path = NautyExecutable("geng").absolute_filename()
            proc = subprocess.run(
                f"{geng_path} -u {n}".split(),
                capture_output=True
            )
            result = re.findall(b">Z ([0-9]+) graphs generated", proc.stderr)
            return Integer(result[0])


class A000108(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=108,
            description="Catalan numbers: C(n) = binomial(2n,n)/(n+1) = (2n)!/(n!(n+1)!).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.combinat.combinat import catalan_number
        return catalan_number(n)


class A000110(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=110,
            description="Bell or exponential numbers: number of ways to partition a set of n labeled elements.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.combinat.combinat import bell_number
        return bell_number(n)


class A000111(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=111,
            description="Euler or up/down numbers: e.g.f. sec(x) + tan(x). Also for n >= 2, half the number of alternating permutations on n letters (A001250).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List[Integer]:
        # We calculate terms using the recurrence
        # a(0) = a(1) = 1
        # a(n) = a(n) = Sum_{i = 0..n-2} binomial(n-2,i)*a(i)*a(n-1-i)
        from sage.functions.other import binomial as C
        a = [1, 1]
        while len(a) < n:
            n_ = len(a)
            a.append(sum([C(n_-2, i)*a[i]*a[n_-1-i] for i in range(n_-1)]))
        return [Integer(i) for i in a[:n]]


class A000112(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=112,
            description="Number of partially ordered sets (\"posets\") with n unlabeled elements.",
            all_at_once=False
        )
        self.cached = cached
        self.SEQ_TERMS = [1, 1, 2, 5, 16, 63, 318, 2045, 16999, 183231, 2567284,
                          46749427, 1104891746, 33823827452, 1338193159771,
                          68275077901156, 4483130665195087]

    def _eval(self, n: Integer) -> Integer:
        if self.cached:
            return Integer(self.SEQ_TERMS[n])
        else:
            from sage.combinat.posets.posets import FinitePosets_n
            return len(FinitePosets_n(n))


class A000120(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=120,
            description="1's-counting sequence: number of 1's in binary expansion of n (or the binary weight of n).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(n).popcount()


class A000123(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=123,
            description="Number of binary partitions: number of partitions of 2n into powers of 2.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Using the recurrence
        # a(n) = a(n-1) + a(n//2)
        a = [1]
        for i in range(1, n):
            a.append(a[i-1]+a[i//2])
        return [Integer(i) for i in a]


class A000124(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=124,
            description="Central polygonal numbers (the Lazy Caterer's sequence): n(n+1)/2 + 1; or, maximal number of pieces formed when slicing a pancake with n cuts.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer((n*(n+1))//2 + 1)


class A000129(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=129,
            description="Pell numbers: a(0) = 0, a(1) = 1; for n > 1, a(n) = 2*a(n-1) + a(n-2).",
            all_at_once=True,
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [0, 1]
        for i in range(2, n):
            a.append(2*a[i-1]+a[i-2])
        return [Integer(i) for i in a[:n]]


class A000140(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=140,
            description="Kendall-Mann numbers: the most common number of inversions in a permutation on n letters is floor(n*(n-1)/4); a(n) is the number of permutations with this many inversions.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        # Using formula due to David W. Wilson
        # Largest coefficient of (1)(x+1)(x^2+x+1)...(x^(n-1) + ... + x + 1)
        from sage.all import PolynomialRing, ZZ
        R = PolynomialRing(ZZ, names=('x',))
        (x,) = R.gens()
        return max(arith.prod([sum([x**j for j in range(i+1)]) for i in range(n)]).coefficients())


class A000142(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=142,
            description="Factorial numbers: n! = 1*2*3*4*...*n (order of symmetric group S_n, number of permutations of n letters).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.functions.other import factorial
        return factorial(n)


class A000161(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=161,
            description="Number of partitions of n into 2 squares.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        # See comment by Ant King
        def f(n): return sum([
            1 if d % 4 == 1 else -1 if d % 4 == 3 else 0
            for d in arith.divisors(n)
        ])
        def delta(n): return 1 if Integer(n).is_square() else 0
        # a(n)=1/2 (f(n)+delta(n)+delta(1/2 n))
        return 1 if n == 0 else (f(n) + delta(n) + delta(n // 2))//2


class A000166(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=166,
            description="Subfactorial or rencontres numbers, or derangements: number of permutations of n elements with no fixed points.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Using recurrence
        # a(n) = n*a(n-1) + (-1)^n
        a = [1]
        for i in range(1, n):
            a.append(i*a[-1] + pow(-1, i))
        return [Integer(i) for i in a[:n]]


class A000169(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=169,
            description="Number of labeled rooted trees with n nodes: n^(n-1).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(pow(n, n-1))


class A000182(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=182,
            description="Tangent (or \"Zag\") numbers: e.g.f. tan(x), also (up to signs) e.g.f. tanh(x).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Using implementation from Peter Luschny
        T = [0, 1] + [0] * (n-1)
        for k in range(2, n+1):
            T[k] = (k-1)*T[k-1]
        for k in range(2, n+1):
            for j in range(k, n+1):
                T[j] = (j-k)*T[j-1]+(j-k+2)*T[j]
        return [Integer(i) for i in T[1:]]


class A000203(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=203,
            description="a(n) = sigma(n), the sum of the divisors of n. Also called sigma_1(n).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return arith.sigma(n, 1)


class A000204(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=204,
            description="Lucas numbers (beginning with 1): L(n) = L(n-1) + L(n-2) with L(1) = 1, L(2) = 3.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = [1, 3]
        while len(seq) < n:
            seq.append(seq[-1] + seq[-2])
        return [Integer(i) for i in seq[:n]]


class A000217(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=217,
            description="Triangular numbers: a(n) = binomial(n+1,2) = n*(n+1)/2 = 0 + 1 + 2 + ... + n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer((n*(n+1))//2)


class A000219(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=219,
            description="Number of planar partitions (or plane partitions) of n.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Using the recurrence
        # a(n) = (1/n) * Sum_{k=1..n} a(n-k)*sigma_2(k), n > 0, a(0)=1
        # From Vladeta Jovovic
        a = [1]
        for n_ in range(1, n+1):
            an = sum([a[n_-k]*arith.sigma(k, 2) for k in range(1, n_+1)])//n_
            a.append(an)
        return [Integer(i) for i in a[:n]]


class A000225(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=225,
            description="a(n) = 2^n - 1. (Sometimes called Mersenne numbers, although that name is usually reserved for A001348.)",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer((1 << n)-1)


class A000244(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=225,
            description="Powers of 3: a(n) = 3^n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(3**n)


class A000262(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=262,
            description="Number of \"sets of lists\": number of partitions of {1,...,n} into any number of lists, where a list means an ordered subset.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Using fact that sequence is D-finite
        a = [1, 1]
        for k in range(2, n+1):
            a.append((2*k-1)*a[k-1] - (k-1)*(k-2)*a[k-2])
        return [Integer(i) for i in a[:n]]


class A000272(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=272,
            description="Number of trees on n labeled nodes: n^(n-2) with a(0)=1.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(1) if n == 0 else Integer(pow(n, n-2))


class A000273(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=273,
            description="Number of unlabeled simple digraphs with n nodes.",
            all_at_once=False,
        )
        self.cached = cached
        self.SEQ_TERMS = [1, 1, 3, 16, 218, 9608, 1540944, 882033440,
                          1793359192848, 13027956824399552,
                          341260431952972580352, 32522909385055886111197440,
                          11366745430825400574433894004224,
                          14669085692712929869037096075316220928,
                          70315656615234999521385506555979904091217920]

    def _eval(self, n: Integer) -> Integer:
        if self.cached:
            return Integer(self.SEQ_TERMS[n])
        else:
            if n == 0:
                return 1
            # Use nauty to count digraphs
            # This is impractical on my machine for any n > 6

            # The equivalent shell command is
            # nauty-geng {n} | nauty-directg -u

            from sage.features.nauty import NautyExecutable
            from subprocess import Popen, PIPE
            import re

            geng_path = NautyExecutable("geng").absolute_filename()
            directg_path = NautyExecutable("directg").absolute_filename()

            geng_proc = Popen(
                f"{geng_path} {n}".split(),
                stdout=PIPE,
                stderr=PIPE
            )
            directg_proc = Popen(
                f"{directg_path} -u".split(),
                stdin=geng_proc.stdout,
                stdout=PIPE,
                stderr=PIPE
            )
            output = [i for i in directg_proc.communicate() if i]
            res = sum([[Integer(i) for i in re.findall(
                b"([0-9]+) digraphs generated", s)] for s in output], [])
            return res[0]


class A000290(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=290,
            description="The squares: a(n) = n^2.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(n*n)


class A000292(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=292,
            description="Tetrahedral (or triangular pyramidal) numbers: a(n) = C(n+2,3) = n*(n+1)*(n+2)/6.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer((n*(n+1)*(n+2))//6)


class A000302(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=302,
            description="Powers of 4: a(n) = 4^n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(pow(4, n))


class A000311(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=311,
            description="Schroeder's fourth problem; also series-reduced rooted trees with n labeled leaves; also number of total partitions of n.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # From given recurrence
        from sage.functions.other import binomial as C
        a = [0, 1, 1]
        for k in range(2, n+1):
            a.append((k+2)*a[k] + 2*sum([C(k, j)*a[j]*a[k-j+1]
                     for j in range(2, k)]))
        return [Integer(i) for i in a[:n]]


class A000312(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=312,
            description="a(n) = n^n; number of labeled mappings from n points to themselves (endofunctions).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(pow(n, n))


class A000326(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=326,
            description="Pentagonal numbers: a(n) = n*(3*n-1)/2.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer((n*(3*n-1)//2))


class A000330(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=330,
            description="Square pyramidal numbers: a(n) = 0^2 + 1^2 + 2^2 + ... + n^2 = n*(n+1)*(2*n+1)/6.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer((n*(n+1)*(2*n+1))//6)


class A000364(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=364,
            description="Euler (or secant or \"Zig\") numbers: e.g.f. (even powers only) sec(x) = 1/cos(x).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.combinat.combinat import euler_number
        return Integer(abs(euler_number(2*n)))


class A000396(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=396,
            description="Perfect numbers k: k is equal to the sum of the proper divisors of k.",
            all_at_once=True
        )
        self.cached = cached

    def _eval_up_to_n(self, n: Integer) -> List:
        if self.cached:
            return [Integer(pow(2, p-1)*(pow(2, p)-1)) for p in A000043()._eval_up_to_n(n)]
        else:
            # Note that the given terms are the only known perfect numbers
            # so the given code below is only for completeness
            seq = []
            for k in count(1):
                if arith.sigma(k, 1) == 2*k:
                    seq.append(k)
                if len(seq) == n:
                    return [Integer(i) for i in seq]


class A000521(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=521,
            description="Coefficients of modular function j as power series in q = e^(2 Pi i t). Another name is the elliptic modular invariant J(tau).",
            all_at_once=True,
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        from sage.modular.modform.j_invariant import j_invariant_qexp
        coeffs = j_invariant_qexp(n).coefficients()
        return [Integer(i) for i in coeffs[:n]]


class A000578(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=578,
            description="The cubes: a(n) = n^3.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(pow(n, 3))


class A000583(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=583,
            description="Fourth powers: a(n) = n^4.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(pow(n, 4))


class A000593(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=593,
            description="Sum of odd divisors of n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(sum([d for d in arith.divisors(n) if d % 2 == 1]))


class A000594(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=594,
            description="Ramanujan's tau function (or Ramanujan numbers, or tau numbers).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        from sage.modular.modform.vm_basis import delta_qexp
        return [Integer(i) for i in list(delta_qexp(n+1))[1:]]


class A000602(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=602,
            description="Number of n-node unrooted quartic trees; number of n-carbon alkanes C(n)H(2n+2) ignoring stereoisomers.",
            all_at_once=False
        )
        self.cached = cached
        self.SEQ_TERMS = [1, 1, 1, 1, 2, 3, 5, 9, 18, 35, 75, 159, 355, 802, 1858, 4347,
                          10359, 24894, 60523, 148284, 366319, 910726, 2278658,
                          5731580, 14490245, 36797588, 93839412, 240215803,
                          617105614, 1590507121, 4111846763, 10660307791,
                          27711253769]

    def _eval(self, n: Integer) -> Integer:
        if self.cached:
            return Integer(self.SEQ_TERMS[n])
        else:
            if n == 0:
                return 1
            # Use the nauty-gentreeg to enumerate trees
            from sage.features.nauty import NautyExecutable
            import subprocess
            import re

            geng_path = NautyExecutable("gentreeg").absolute_filename()
            proc = subprocess.run(
                f"{geng_path} {n} -u -D4".split(),
                capture_output=True
            )
            result = re.findall(b">Z ([0-9]+) trees generated", proc.stderr)
            return Integer(result[0])


class A000670(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=670,
            description="Fubini numbers: number of preferential arrangements of n labeled elements; or number of weak orders on n labeled elements; or number of ordered partitions of [n].",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Using recurrence
        # a(n) = Sum_{k=1..n} binomial(n, k)*a(n-k), a(0) = 1.
        # a(k) = Sum_{j=1..k} binomial(k, j)*a(k-j), a(0) = 1.
        from sage.functions.other import binomial as C
        a = [1]
        for k in range(1, n+1):
            a.append(sum([a[k-j]*C(k, j) for j in range(1, k+1)]))
        return [Integer(i) for i in a[:n]]


class A000688(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=688,
            description="Number of Abelian groups of order n; number of factorizations of n into prime powers.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.combinat.partition import Partitions
        return arith.prod(Partitions(i).cardinality() for i in dict(arith.factor(n)).values())


class A000720(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=720,
            description="pi(n), the number of primes <= n. Sometimes called PrimePi(n) to distinguish it from the number 3.14159...",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.functions.prime_pi import prime_pi
        return Integer(prime_pi(n))


class A000793(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=793,
            description="Landau's function g(n): largest order of permutation of n elements. Equivalently, largest LCM of partitions of n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.arith.functions import lcm
        from sage.combinat.partition import Partitions
        return max(lcm(l) for l in Partitions(n))


class A000796(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=796,
            description="Decimal expansion of Pi (or digits of Pi).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        from sage.all import pi
        return [Integer(i) for i in str(pi.n(digits=n+10)) if i in "0123456789"][:n]


class A000959(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=959,
            description="Lucky numbers.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        if n > 200000:
            raise ValueError("n too large")
        # Code from Robert FERREOL from OEIS page
        L = list(range(1, max(n*n, 10**6), 2))
        j = 1
        while j <= len(L) - 1 and L[j] <= len(L):
            del L[L[j]-1::L[j]]
            j += 1
        return [Integer(i) for i in L[:n]]


class A000961(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=961,
            description="Powers of primes. Alternatively, 1 and the prime powers (p^k, p prime, k >= 1).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count(1):
            if i == 1 or Integer(i).is_prime_power():
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq][:n]


class A000984(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=984,
            description="Central binomial coefficients: binomial(2*n,n) = (2*n)!/(n!)^2.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.functions.other import binomial as C
        return Integer(C(2*n, n))


class A001003(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1003,
            description="Schroeder's second problem (generalized parentheses); also called super-Catalan numbers or little Schroeder numbers.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Using recurrence
        # (n+1) * a(n) = (6*n-3) * a(n-1) - (n-2) * a(n-2) if n>1. a(0) = a(1) = 1.
        a = [1, 1]
        for k in range(2, n+1):
            a.append(((6*k-3)*a[k-1] - (k-2)*a[k-2]) // (k+1))
        return [Integer(i) for i in a[:n]]

class A001006(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1006,
            description="Motzkin numbers: number of ways of drawing any number of nonintersecting chords joining n (labeled) points on a circle.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        # Using recurrence (n+2)*a(n) = (2*n+1)*a(n-1) + (3*n-3)*a(n-2).
        a = [1,1]
        for k in range(2, n+1):
            a.append( ( (2*k+1)*a[k-1] + (3*k-3)*a[k-2] ) // (k+2))
        return [Integer(i) for i in a[:n]]

class A001037(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1037,
            description="Number of degree-n irreducible polynomials over GF(2); number of n-bead necklaces with beads of 2 colors when turning over is not allowed and with primitive period n; number of binary Lyndon words of length n.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(1) if n == 0 else Integer(sum([arith.moebius(n//d)*pow(2,d) for d in arith.divisors(n)])//n)
    
class A001045(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1045,
            description="Jacobsthal sequence (or Jacobsthal numbers): a(n) = a(n-1) + 2*a(n-2), with a(0) = 0, a(1) = 1; also a(n) = nearest integer to 2^n/3.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer( (pow(2,n)-(pow(-1,n)))//3 )

class A001055(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1055,
            description="The multiplicative partition function: number of ways of factoring n with all factors greater than 1 (a(1) = 1 by convention).",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        # Translation of Python code provided on OEIS page
        def T(n,m):
            if arith.is_prime(n): return 1 if n<=m else 0
            s = sum([T(n//d, d) for d in arith.divisors(n) if 1 < d <= min(m,n-1)])
            return s + (1 if n<=m else 0)
        return Integer(T(n,n))
    
class A001057(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1057,
            description="Canonical enumeration of integers: interleaved positive and negative integers with zero prepended.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer( ( 1-(2*n+1)*(pow(-1,n)) ) // 4)

class A001065(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1065,
            description="Sum of proper divisors (or aliquot parts) of n: sum of divisors of n that are less than n.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.sigma(n, 1) - n)
    
class A001097(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1097,
            description="Twin primes.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        seq = set()
        p, q = 2, 3
        while len(seq) < n:
            if q - p == 2: 
                seq |= {p,q}
            p, q = q, arith.next_prime(q)
        return [Integer(i) for i in sorted(seq)]

class A001113(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1113,
            description="Decimal expansion of e.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        from sage.all import e
        return [Integer(i) for i in str(e.n(digits=n+10)) if i in "0123456789"][:n]
    
class A001147(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1147,
            description="Double factorial of odd numbers: a(n) = (2*n-1)!! = 1*3*5*...*(2*n-1).",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.rising_factorial(n+1,n)//pow(2, n))
    
class A001157(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1157,
            description="a(n) = sigma_2(n): sum of squares of divisors of n.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return arith.sigma(n, 2)

class A001190(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1190,
            description="Wedderburn-Etherington numbers: unlabeled binary rooted trees (every node has outdegree 0 or 2) with n endpoints (and 2n-1 nodes in all).",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        # Use recurrence given
        # Use functools.cache instead of array for recurrence due to more complex nature of recurrence
        @cache
        def a(n):
            if n == 0: return 0
            if n == 1: return 1
            m = n//2+(1 if n%2==1 else 0)
            if n % 2 == 1: # n = 2*m-1
                return sum([a(i)*a(2*m-1-i) for i in range(1,m)])
            else: # n = 2*m
                return sum([a(i)*a(2*m-i) for i in range(1,m)]) + (a(m)*(a(m)+1))//2
        return [Integer(a(i)) for i in range(n)]

class A001221(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1221,
            description="Number of distinct primes dividing n (also called omega(n)).",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(len(arith.prime_divisors(n)))
    
class A001222(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1222,
            description="Number of prime divisors of n counted with multiplicity (also called big omega of n, bigomega(n) or Omega(n)).",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(sum(e for _, e in dict(arith.factor(n)).items()))

class A001227(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1227,
            description="Number of odd divisors of n.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(sum([1 for d in arith.divisors(n) if d % 2 == 1]))
    
class A001285(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1285,
            description="Thue-Morse sequence: let A_k denote the first 2^k terms; then A_0 = 1 and for k >= 0, A_{k+1} = A_k B_k, where B_k is obtained from A_k by interchanging 1's and 2's.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(Integer(n).popcount() % 2 + 1)