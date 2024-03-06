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
Core sequences to implement:
 - A000014: Number of series-reduced trees with n nodes.
 - A000105: Number of free polyominoes (or square animals) with n cells.
 - A000109: Number of simplicial polyhedra with n vertices; simple planar graphs with n vertices and 3n-6 edges; maximal simple planar graphs with n vertices; planar triangulations with n vertices; triangulations of the sphere with n vertices; 3-connected cubic planar graphs on 2n-4 vertices.
 - A000609: Number of threshold functions of n or fewer variables.
 - A000612: Number of P-equivalence classes of switching functions of n or fewer variables, divided by 2.
 - A000798: Number of different quasi-orders (or topologies, or transitive digraphs) with n labeled elements.
 - A001034: Orders of noncyclic simple groups (without repetition).
 - A002572: Number of partitions of 1 into n powers of 1/2; or (according to one definition of "binary") the number of binary rooted trees.
 - A005470: Number of unlabeled planar simple graphs with n nodes.
 - A005588: Number of free binary trees admitting height n.
  - paper with recurrence: https://oeis.org/A005588/a005588.pdf
 - A055512: Lattices with n labeled elements.
 - A104725: Number of complementing systems of subsets of {0, 1, ..., n-1}. 
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
        a = [1, 1]
        while len(a) < n:
            n_ = len(a)
            a.append(sum([arith.binomial(n_-2, i)*a[i]*a[n_-1-i]
                     for i in range(n_-1)]))
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
            seq_number=244,
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
        a = [0, 1, 1]
        for k in range(2, n+1):
            a.append((k+2)*a[k] + 2*sum([arith.binomial(k, j)*a[j]*a[k-j+1]
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
        a = [1]
        for k in range(1, n+1):
            a.append(sum([a[k-j]*arith.binomial(k, j) for j in range(1, k+1)]))
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
        return Integer(arith.binomial(2*n, n))


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
        a = [1, 1]
        for k in range(2, n+1):
            a.append(((2*k+1)*a[k-1] + (3*k-3)*a[k-2]) // (k+2))
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
        return Integer(1) if n == 0 else Integer(sum([arith.moebius(n//d)*pow(2, d) for d in arith.divisors(n)])//n)


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
        return Integer((pow(2, n)-(pow(-1, n)))//3)


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
        def T(n, m):
            if arith.is_prime(n):
                return 1 if n <= m else 0
            s = sum([T(n//d, d)
                    for d in arith.divisors(n) if 1 < d <= min(m, n-1)])
            return s + (1 if n <= m else 0)
        return Integer(T(n, n))


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
        return Integer((1-(2*n+1)*(pow(-1, n))) // 4)


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
                seq |= {p, q}
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
        return Integer(2*n-1).multifactorial(2)


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
            if n == 0:
                return 0
            if n == 1:
                return 1
            m = n//2+(1 if n % 2 == 1 else 0)
            if n % 2 == 1:  # n = 2*m-1
                return sum([a(i)*a(2*m-1-i) for i in range(1, m)])
            else:  # n = 2*m
                return sum([a(i)*a(2*m-i) for i in range(1, m)]) + (a(m)*(a(m)+1))//2
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


class A001333(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1333,
            description="Numerators of continued fraction convergents to sqrt(2).",
            all_at_once=True,
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [1, 1]
        while len(a) < n:
            a.append(2*a[-1] + a[-2])
        return [Integer(i) for i in a][:n]


class A001349(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1349,
            description="Number of connected graphs with n nodes.",
            all_at_once=False
        )
        self.cached = cached
        self.SEQ_TERMS = [1, 1, 1, 2, 6, 21, 112, 853, 11117, 261080, 11716571,
                          1006700565, 164059830476, 50335907869219,
                          29003487462848061, 31397381142761241960,
                          63969560113225176176277,
                          245871831682084026519528568,
                          1787331725248899088890200576580,
                          24636021429399867655322650759681644]

    def _eval(self, n: Integer) -> Integer:
        if self.cached:
            return Integer(self.SEQ_TERMS[n])
        else:
            # Use nauty-geng to enumerate all connected graphs on n nodes
            if n == 0:
                return 1
            from sage.features.nauty import NautyExecutable
            import subprocess
            import re

            geng_path = NautyExecutable("geng").absolute_filename()
            proc = subprocess.run(
                f"{geng_path} {n} -c -u".split(),
                capture_output=True
            )
            result = re.findall(b">Z ([0-9]+) graphs generated", proc.stderr)
            return Integer(result[0])


class A001358(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1358,
            description="Semiprimes (or biprimes): products of two primes.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count(2):
            if sum(dict(arith.factor(i)).values()) == 2:
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]


class A001405(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1405,
            description="a(n) = binomial(n, floor(n/2)).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.binomial(n, n//2))


class A001462(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1462,
            description="Golomb's sequence: a(n) is the number of times n occurs, starting with a(1) = 1.",
            all_at_once=True,
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [None, 1, 2, 2]
        k = 3
        while len(a) < n+1:
            a += [k] * a[k]
            k += 1
        return [Integer(i) for i in a[1:n+1]]


class A001477(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1477,
            description="The nonnegative integers.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(n)


class A001478(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1478,
            description="The negative integers.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(-n)


class A001481(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1481,
            description="Numbers that are the sum of 2 squares.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        from sage.rings.sum_of_squares import is_sum_of_two_squares_pyx
        seq = []
        for i in count(0):
            if is_sum_of_two_squares_pyx(i):
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]


class A001489(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1489,
            description="a(n) = -n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(-n)


class A001511(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1511,
            description='The ruler function: exponent of the highest power of 2 dividing 2n. Equivalently,  the 2-adic valuation of 2n.',
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.valuation(2*n, 2))


class A001519(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1519,
            description="a(n) = 3*a(n-1) - a(n-2) for n >= 2, with a(0) = a(1) = 1.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [1, 1]
        while len(a) < n:
            a.append(3*a[-1] - a[-2])
        return [Integer(i) for i in a]


class A001615(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1615,
            description="Dedekind psi function: n * Product_{p|n, p prime} (1 + 1/p).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(n*arith.prod([(p+1)/p for p in arith.prime_divisors(n)]))


class A001699(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1699,
            description="Number of binary trees of height n; or products (ways to insert parentheses) of height n when multiplication is non-commutative and non-associative.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        from sage.misc.functional import isqrt
        a = [1, 1]
        while len(a) < n:
            a.append(a[-1]**2 + a[-1] + a[-1]*isqrt(4*a[-1]-3))
        return [Integer(i) for i in a]


class A001700(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1700,
            description="a(n) = binomial(2*n+1, n+1): number of ways to put n+1 indistinguishable balls into n+1 distinguishable boxes = number of (n+1)-st degree monomials in n+1 variables = number of monotone maps from 1..n+1 to 1..n+1.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.binomial(2*n+1, n+1))


class A001764(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1764,
            description="a(n) = binomial(3*n,n)/(2*n+1) (enumerates ternary trees and also noncrossing trees).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.binomial(3*n, n) // (2*n+1))


class A001906(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=1906,
            description="F(2n) = bisection of Fibonacci sequence: a(n) = 3*a(n-1) - a(n-2).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [0, 1]
        while len(a) < n:
            a.append(3*a[-1]-a[-2])
        return [Integer(i) for i in a]


class A001969(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=1969,
            description="Evil numbers: nonnegative integers with an even number of 1's in their binary expansion.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count():
            if Integer(i).popcount() % 2 == 0:
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]


class A002033(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=2033,
            description="Number of perfect partitions of n.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [1, 1]
        for i in range(2, n):
            a.append(sum([a[j-1] for j in arith.divisors(i+1) if j <= i]))
        return [Integer(i) for i in a]


class A002083(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=2083,
            description="Narayana-Zidek-Capell numbers: a(n) = 1 for n <= 2. Otherwise a(2n) = 2a(2n-1), a(2n+1) = 2a(2n) - a(n).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [None, 1, 1]
        for k in range(3, n+1):
            a.append(2*a[k-1] - (a[k//2] if k % 2 == 1 else 0))
        return [Integer(i) for i in a[1:]]


class A002106(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=2106,
            description="Number of transitive permutation groups of degree n.",
            all_at_once=False
        )
        self.cached = cached
        self.SEQ_TERMS = [None, 1, 1, 2, 5, 5, 16, 7, 50, 34, 45, 8, 301, 9, 63, 104, 1954, 10,
                          983, 8, 1117, 164, 59, 7, 25000, 211, 96, 2392, 1854, 8, 5712,
                          12, 2801324, 162, 115, 407, 121279, 11, 76, 306, 315842, 10,
                          9491, 10, 2113, 10923, 56, 6]

    def _eval(self, n: Integer) -> Integer:
        if self.cached:
            return Integer(self.SEQ_TERMS[n])
        else:
            from sage.libs.gap.libgap import libgap
            return Integer(len(libgap.AllTransitiveGroups(libgap.NrMovedPoints, n)))


class A002110(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=2110,
            description="Primorial numbers (first definition): product of first n primes. Sometimes written prime(n)#.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.prod(arith.primes_first_n(n)))


class A002113(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=2113,
            description="Palindromes in base 10.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count():
            if str(i)[::-1] == str(i):
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]


class A002275(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=2275,
            description="Repunits: (10^n - 1)/9. Often denoted by R_n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer((pow(10, n)-1)//9)


class A002322(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=2322,
            description="Reduced totient function psi(n): least k such that x^k == 1 (mod n) for all x prime to n; also known as the Carmichael lambda function (exponent of unit group mod n); also called the universal exponent of n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.carmichael_lambda(n))


class A002378(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=2378,
            description="Oblong (or promic, pronic, or heteromecic) numbers: a(n) = n*(n+1).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(n*(n+1))


class A002426(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=2426,
            description="Central trinomial coefficients: largest coefficient of (1 + x + x^2)^n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(sum([arith.binomial(n, k)*arith.binomial(k, n-k) for k in range(n+1)]))


class A002487(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=2487,
            description="Stern's diatomic series (or Stern-Brocot sequence): a(0) = 0, a(1) = 1; for n > 0: a(2*n) = a(n), a(2*n+1) = a(n) + a(n+1).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [0, 1]
        for k in range(2, n):
            a.append(a[k//2]+(a[k//2+1] if k % 2 == 1 else 0))
        return [Integer(i) for i in a]


class A002530(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=2530,
            description="a(n) = 4*a(n-2) - a(n-4) for n > 1, a(n) = n for n = 0, 1.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [0, 1, 1, 3]
        while len(a) < n:
            a.append(4*a[-2]-a[-4])
        return [Integer(i) for i in a]


class A002531(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=2531,
            description="a(2*n) = a(2*n-1) + a(2*n-2), a(2*n+1) = 2*a(2*n) + a(2*n-1); a(0) = a(1) = 1.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = [1, 1]
        for i in range(2, n):
            if i % 2 == 0:
                seq.append(seq[-1]+seq[-2])
            else:
                seq.append(2*seq[-1]+seq[-2])
        return [Integer(i) for i in seq]


class A002620(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=2620,
            description="Quarter-squares: a(n) = floor(n/2)*ceiling(n/2). Equivalently, a(n) = floor(n^2/4).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer((n*n)//4)


class A002654(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=2654,
            description="Number of ways of writing n as a sum of at most two nonzero squares, where order matters; also (number of divisors of n of form 4m+1) - (number of divisors of form 4m+3).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(sum([
            1 if d % 4 == 1 else -1 if d % 4 == 3 else 0
            for d in arith.divisors(n)
        ]))


class A002658(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=2658,
            description="a(0) = a(1) = 1; for n > 0, a(n+1) = a(n)*(a(0) + ... + a(n-1)) + a(n)*(a(n) + 1)/2.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        a = [1, 1]
        for k in range(1, n-1):
            a.append(a[k]*sum([a[i] for i in range(k)]) + (a[k]*(a[k]+1))//2)
        return [Integer(i) for i in a]


class A002808(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=2808,
            description="The composite numbers: numbers n of the form x*y for x > 1 and y > 1.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count(2):
            if not arith.is_prime(i):
                seq.append(Integer(i))
            if len(seq) == n:
                return seq


class A003094(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=3094,
            description="Number of unlabeled connected planar simple graphs with n nodes.",
            all_at_once=False
        )
        self.cached = cached
        self.SEQ_TERMS = [1, 1, 1, 2, 6, 20, 99, 646, 5974,
                          71885, 1052805, 17449299, 313372298, 5942258308]

    def _eval(self, n: Integer) -> Integer:
        if self.cached:
            return Integer(self.SEQ_TERMS[n])
        else:
            # Using provided nauty program
            #  geng -c $n | planarg -q | countg -q
            if n == 0:
                return 1
            from sage.features.nauty import NautyExecutable
            from subprocess import Popen, PIPE
            import re
            # NautyExecutable("geng").absolute_filename()
            geng_proc = Popen(
                f"{NautyExecutable('geng').absolute_filename()} -c {n}".split(),
                stdout=PIPE,
                stderr=PIPE
            )
            planarg_proc = Popen(
                f"{NautyExecutable('planarg').absolute_filename()} -q".split(),
                stdin=geng_proc.stdout,
                stdout=PIPE,
                stderr=PIPE
            )
            countg_proc = Popen(
                f"{NautyExecutable('countg').absolute_filename()} -q".split(),
                stdin=planarg_proc.stdout,
                stdout=PIPE,
                stderr=PIPE
            )
            output = countg_proc.communicate()[0]
            return Integer(re.findall(b"([0-9]+) graphs", output)[0])


class A003136(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=3136,
            description="Loeschian numbers: numbers of the form x^2 + xy + y^2; norms of vectors in A2 lattice.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        def _is_ok(n):
            for p, e in dict(arith.factor(n)).items():
                if p % 3 == 2 and e % 2 == 1:
                    return False
            return True
        seq = [0]
        for i in count(1):
            if _is_ok(i):
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]


class A003418(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=3418,
            description="Least common multiple (or LCM) of {1, 2, ..., n} for n >= 1, a(0) = 1.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        from sage.arith.functions import lcm
        return Integer(1) if n == 0 else Integer(lcm(range(1, n+1)))


class A003484(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=3484,
            description="Radon function, also called Hurwitz-Radon numbers.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(8*(arith.valuation(n, 2)//4)+pow(2, arith.valuation(n, 2) % 4))


class A004011(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=4011,
            description="Theta series of D_4 lattice; Fourier coefficients of Eisenstein series E_{gamma,2}.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(1) if n == 0 else 24*sum([d for d in arith.divisors(n) if d % 2 == 1])


class A004018(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=4018,
            description="Theta series of square lattice (or number of ways of writing n as a sum of 2 squares). Often denoted by r(n) or r_2(n).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(1) if n == 0 else Integer(
            4*arith.prod([
                (1+e) if (p % 4 == 1) else
                (0 if e % 2 == 1 else 1) if (p % 4 == 3) else 1
                for p, e in dict(arith.factor(n)).items()
            ])
        )


class A004526(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=4526,
            description="Nonnegative integers repeated, floor(n/2).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(n//2)


class A005036(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=5036,
            description="Number of nonequivalent dissections of a polygon into n quadrilaterals by nonintersecting diagonals up to rotation and reflection.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        # Implementation from this paper: https://arxiv.org/pdf/1807.11602.pdf
        def v(m):
            # Equation 2.1
            return arith.binomial(3*m, m)/(2*m+1)

        def s(m):
            # Theorem 2.8
            k = m // 2
            if m % 2 == 0:
                # k = 2*m
                return arith.binomial(3*k, k) / (2*k+1)
            else:
                # k = 2*m + 1
                return arith.binomial(3*k+1, k) / (k+1)

        def a(m):
            assert m % 2 == 0
            n = m//2
            # Theorem 4.2
            if n % 2 == 1:
                return (v(n-1) + 3*n*s(n-1)) / (4*n)
            elif n % 4 == 0:
                return (v(n-1) + 5*n*s(n-1)/2) / (4*n)
            elif n % 4 == 2:
                return (v(n-1) + 5*n*s(n-1)/2 + n*s((n-2)/2)) / (4*n)
        return Integer(a(2*(n+1)))


class A005100(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=5100,
            description="Deficient numbers: numbers k such that sigma(k) < 2k.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count(1):
            if sum(arith.divisors(i)) < 2*i:
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]


class A005101(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=5101,
            description="Abundant numbers (sum of divisors of m exceeds 2m).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count(1):
            if sum(arith.divisors(i)) > 2*i:
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]


class A005117(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=5117,
            description="Squarefree numbers: numbers that are not divisible by a square greater than 1.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count(1):
            if arith.is_squarefree(i):
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]


class A005130(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=5130,
            description="Robbins numbers: a(n) = Product_{k=0..n-1} (3k+1)!/(n+k)!; also the number of descending plane partitions whose parts do not exceed n; also the number of n X n alternating sign matrices (ASM's).",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return arith.prod([arith.factorial(3*k+1) / arith.factorial(n+k) for k in range(n)])


class A005230(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=5230,
            description="Stern's sequence: a(1) = 1, a(n+1) is the sum of the m preceding terms, where m*(m-1)/2 < n <= m*(m+1)/2 or equivalently m = ceiling((sqrt(8*n+1)-1)/2) = A002024(n).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        from math import isqrt
        seq = [None, 1]
        for k in range(2, n+1):
            m = (isqrt(8*(k-1))+1)//2
            seq.append(sum(seq[-m:]))
        return [Integer(i) for i in seq[1:]]


class A005408(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=5408,
            description="The odd numbers: a(n) = 2*n + 1.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(2*n+1)


class A005811(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=5811,
            description="Number of runs in binary expansion of n (n>0); number of 1's in Gray code for n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(n ^ (n//2)).popcount()


class A005843(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=5843,
            description="The nonnegative even numbers: a(n) = 2n.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(n*2)


class A006318(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=6318,
            description="Large Schröder numbers (or large Schroeder numbers, or big Schroeder numbers).",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Using recurrence (n-2)*a(n-2) - 3*(2*n-1)*a(n-1) + (n+1)*a(n) = 0
        # => (n-2)*a(n-2) - 3*(2*n-1)*a(n-1) = -(n+1)*a(n)
        seq = [1, 2]
        for k in range(2, n):
            seq.append(((k-2)*seq[k-2] - 3*(2*k-1)*seq[k-1]) / (-(k+1)))
        return [Integer(i) for i in seq]


class A006530(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=6530,
            description="Gpf(n): greatest prime dividing n, for n >= 2; a(1)=1.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(1) if n == 1 else max(arith.prime_divisors(n))


class A006882(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=6882,
            description="Double factorials n!!: a(n) = n*a(n-2) for n > 1, a(0) = a(1) = 1.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return Integer(n).multifactorial(2)


class A006894(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=6894,
            description="Number of planted 3-trees of height < n.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = [1]
        while len(seq) < n:
            seq.append((seq[-1] * (seq[-1]+1))//2 + 1)
        return [Integer(i) for i in seq]


class A006966(OEISSequence):
    def __init__(self, cached=True):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=6966,
            description="Number of lattices on n unlabeled nodes.",
            all_at_once=False
        )
        self.cached = cached
        self.SEQ_TERMS = [1, 1, 1, 1, 2, 5, 15, 53, 222, 1078, 5994, 37622, 262776,
                          2018305, 16873364, 152233518, 1471613387, 15150569446,
                          165269824761, 1901910625578, 23003059864006]

    def _eval(self, n: Integer) -> Integer:
        if self.cached:
            return Integer(self.SEQ_TERMS[n])
        else:
            from sage.combinat.posets.posets import FinitePosets_n
            return Integer(sum([i.is_lattice() for i in Posets(6)]))

class A007318(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=7318,
            description="Pascal's triangle read by rows: C(n,k) = binomial(n,k) = n!/(k!*(n-k)!), 0 <= k <= n.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        def _table_iterator():
            for n_ in count():
                for k_ in range(n_+1):
                    yield arith.binomial(n_, k_)
        return [Integer(i) for i in islice(_table_iterator(), n)]

class A008275(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=8275,
            description="Triangle read by rows of Stirling numbers of first kind, s(n,k), n >= 1, 1 <= k <= n.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        from sage.combinat.combinat import stirling_number1
        def _table_iterator():
            for n_ in count(1):
                for k_ in range(1, n_+1):
                    # Sequence entries are signed
                    # See https://en.wikipedia.org/wiki/Stirling_numbers_of_the_first_kind#Signs
                    yield pow(-1, n_-k_)*stirling_number1(n_, k_)
        return [Integer(i) for i in islice(_table_iterator(), n)]

class A008277(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=8277,
            description="Triangle of Stirling numbers of the second kind, S2(n,k), n >= 1, 1 <= k <= n.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        from sage.combinat.combinat import stirling_number2
        def _table_iterator():
            for n_ in count(1):
                for k_ in range(1, n_+1):
                    yield stirling_number2(n_, k_)
        return [Integer(i) for i in islice(_table_iterator(), n)]

class A008279(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=8279,
            description="Triangle T(n,k) = n!/(n-k)! (0 <= k <= n) read by rows, giving number of permutations of n things k at a time.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        def _table_iterator():
            for n_ in count():
                for k_ in range(n_+1):
                    yield arith.factorial(n_) // arith.factorial(n_-k_)
        return [Integer(i) for i in islice(_table_iterator(), n)]

class A008292(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=8292,
            description="Triangle of Eulerian numbers T(n,k) (n >= 1, 1 <= k <= n) read by rows.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        def _table_iterator():
            from sage.combinat.combinat import eulerian_number
            for n_ in count(1):
                for k_ in range(1, n_+1):
                    yield 1 if n_ == k_ else eulerian_number(n_, k_)
        return [1] + [Integer(i) for i in islice(_table_iterator(), n-1)]

class A008683(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=8683,
            description="Möbius (or Moebius) function mu(n). mu(1) = 1; mu(n) = (-1)^k if n is the product of k different primes; otherwise mu(n) = 0.",
            all_at_once=False
        )

    def _eval(self, n: Integer) -> Integer:
        return arith.moebius(n)

class A010060(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=10060,
            description="Thue-Morse sequence: let A_k denote the first 2^k terms; then A_0 = 0 and for k >= 0, A_{k+1} = A_k B_k, where B_k is obtained from A_k by interchanging 0's and 1's.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(n).popcount() % 2

class A018252(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=18252,
            description="The nonprime numbers: 1 together with the composite numbers, A002808.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = [1]
        for i in count(2):
            if not arith.is_prime(i):
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]
                

class A020639(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=20639,
            description="Lpf(n): least prime dividing n (when n > 1); a(1) = 1. Or, smallest prime factor of n, or smallest prime divisor of n.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(1) if n == 1 else Integer(min(arith.prime_divisors(n)))

class A020652(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=20652,
            description="Numerators in canonical bijection from positive integers to positive rationals.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        # Code from OEIS page
        # Original author: Indranil Ghosh
        s,k = 0,2
        while s < n: 
            s += arith.euler_phi(k)
            k += 1
        s -= arith.euler_phi(k-1)
        j = 1
        while s < n:
            if arith.gcd(j, k-1) == 1: s += 1
            j += 1
        return Integer(j-1)

class A020653(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=20653,
            description="Denominators in a certain bijection from positive integers to positive rationals.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        # Code from OEIS page
        # Original author: Indranil Ghosh
        s,k = 0,2
        while s < n: 
            s += arith.euler_phi(k)
            k += 1
        s -= arith.euler_phi(k-1)
        j = 1
        while s < n:
            if arith.gcd(j, k-1) == 1: s += 1
            j += 1
        return Integer(k-j)

class A025487(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=25487,
            description="Least integer of each prime signature A124832; also products of primorial numbers A002110.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        def _is_ok(k):
            exponents = [ arith.valuation(k, p) for p in arith.primes(max(arith.prime_divisors(k))+1) ]
            return list(sorted(exponents)) == list(reversed(exponents))
        seq = [1]
        for i in count(2):
            if _is_ok(i):
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]

class A027641(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=27641,
            description="Numerator of Bernoulli number B_n.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.bernoulli(n).numerator())
    
class A027642(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=27642,
            description="Denominator of Bernoulli number B_n.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(arith.bernoulli(n).denominator())
    
class A035099(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=-1,
            seq_number=35099,
            description="McKay-Thompson series of class 2B for the Monster group with a(0) = 40.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        # Translation of Mathematica code by Jean-François Alcover
        from sage.all import PowerSeriesRing, ZZ
        R = PowerSeriesRing(ZZ,default_prec=100,sparse=True,names=('x',))
        (x,) = R.gens()

        out = 1/x
        for k in range(1, n+1):
            out *= 1+x**k
        out = pow(out, -24)

        seq = out.coefficients()
        seq[1] = 40
        return [Integer(i) for i in seq[:n]]

class A038566(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=38566,
            description="Numerators in canonical bijection from positive integers to positive rationals <= 1: arrange fractions by increasing denominator then by increasing numerator.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        def gen_seq():
            yield 1
            for i in count(2):
                for j in range(1, i+1):
                    if arith.gcd(i,j) == 1: yield j
        return [Integer(i) for i in islice(gen_seq(), n)]
    
class A038567(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=38567,
            description="Denominators in canonical bijection from positive integers to positive rationals <= 1.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        totient_sums = [0]
        seq = [None for _ in range(n)]
        for k in count(1):
            totient_sums.append(totient_sums[-1] + arith.euler_phi(k))
            for j in range(totient_sums[k-1], min(n, totient_sums[k])):
                seq[j] = k
            if totient_sums[-1] > n: break
        return [Integer(i) for i in seq]
    
class A038568(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=38568,
            description="Numerators in canonical bijection from positive integers to positive rationals.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        # Translation of the PARI/GP program
        # The rest of the programs seem to give wrong values?
        for q in count(1):
            e = arith.euler_phi(q)
            if n+1 < 2*e:
                for p in count(1):
                    if arith.gcd(p, q) == 1:
                        if n <= 0: return Integer([q, p][n])
                        else: n -= 2
            n -= 2*e
    
class A038569(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=38569,
            description="Denominators in a certain bijection from positive integers to positive rationals.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        # Translation of Python code from Indranil Ghosh
        s,k = 1,2
        while s <= n:
            s += 2*arith.euler_phi(k); k += 1
        s -= 2*arith.euler_phi(k-1)
        j = 1
        while s <= n:
            if arith.gcd(j, k-1) == 1: s += 2
            j += 1
        return Integer(k-1) if s > n+1 else Integer(j-1)
    
class A049310(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=49310,
            description="Triangle of coefficients of Chebyshev's S(n,x) := U(n,x/2) polynomials (exponents in increasing order).",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        # {T(n, k) = if( k<0 || k>n || (n + k)%2, 0, (-1)^((n + k)/2 + k) * binomial((n + k)/2, k))}
        T = lambda n, k: 0 if (k < 0 or k > n or (n+k)%2 == 1) else pow(-1, (n+k)//2+k)*arith.binomial((n+k)//2, k)
        def _table_iterator():
            for n_ in count():
                for k_ in range(n_+1):
                    yield T(n_, k_)
        return [Integer(i) for i in islice(_table_iterator(),n)]

class A070939(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=70939,
            description="Length of binary representation of n.",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        return Integer(1) if n == 0 else Integer(n).bit_length()

class A074206(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=74206,
            description="Kalmár's [Kalmar's] problem: number of ordered factorizations of n.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        seq = [0,1] + [None] * (n-2)
        for k in range(2, n):
            seq[k] = sum([seq[d] for d in arith.divisors(k) if d != k])
        return [Integer(i) for i in seq]

class A217831(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=0,
            seq_number=217831,
            description="Euclid's triangle read by rows. T(n, k) = 1 if k is prime to n, otherwise 0.",
            all_at_once=True
        )
    def _eval_up_to_n(self, n: Integer) -> List:
        def _table_iterator():
            for n_ in count():
                for k_ in range(n_+1):
                    yield 1 if arith.gcd(n_, k_) == 1 else 0
        return [Integer(i) for i in islice(_table_iterator(), n)]

class A226898(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=226898,
            description="Hooley's Delta function: maximum number of divisors of n in [u, eu] for all u. (Here e is Euler's number 2.718... = A001113.)",
            all_at_once=False
        )
    def _eval(self, n: Integer) -> Integer:
        from math import exp
        d = arith.divisors(n); m = 1
        for i in range(len(d)):
            t = exp(1)*d[i]
            m = max(sum([int(d[j] < t) for j in range(i, len(d))]), m)
        return Integer(m)

class A246655(OEISSequence):
    def __init__(self):
        OEISSequence.__init__(
            self,
            offset=1,
            seq_number=246655,
            description="Prime powers: numbers of the form p^k where p is a prime and k >= 1.",
            all_at_once=True
        )

    def _eval_up_to_n(self, n: Integer) -> List:
        seq = []
        for i in count(2):
            if arith.is_prime_power(i):
                seq.append(i)
            if len(seq) == n:
                return [Integer(i) for i in seq]