# First current-improvement oracle: representation repair

The first execution used source SHA256
`15d724ee23052491c67082e621e38e12f81f4e3e0da27aa35d883158135e0967`;
its complete stdout SHA256 is
`3db622849fcc028a59afd3ca507ef975ca36628cd58b15fd67f371a22b8435e7`.
Five predicates passed before the spin-row assertion stopped execution.

The matrix predicate compared algebraically identical but structurally
different SymPy forms: a factored `-q*(a-b)` and expanded `-q*a+q*b`.
The repaired spin and couple predicates simplify their matrix difference
before comparing to zero. No equation, convention, target or physical
input changed. The captured repaired execution passes all eight checks
with process status zero. This implementation/representation failure is
not a scientific verdict about the current transformation.
