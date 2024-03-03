# OEIS sequences in Sagemath

Implementations of various [OEIS sequences](https://oeis.org/) in [SageMath](https://www.sagemath.org/)

The current goal is to implement all sequences with the keyword `core` - there are 179 such sequences and currently around half are implemented. Emphasis is placed on having an implementation that computes terms instead of looking up terms. This may eventually serve as a replacement to (and was heavily inspired by) the `sage.combinat.sloane_functions` library ([linked here](https://github.com/sagemath/sage/blob/develop/src/sage/combinat/sloane_functions.py))

## TODO
 - more comprehensive tests (test class fields, test if cached=True results in same terms as cached=False)
 - comprehensive profiling -> more efficient implementations?
 - more sequences?
 - clean up interface and add support for list slicing
