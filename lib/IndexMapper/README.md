## Index Mapping / Demapping

To take full advantage of IM, we must be able to generate all possible permutations of the available signal constellations, which leads to a factorial number of combinations. If we want to implement this with an ordinary look-up table, the mapping procedure becomes impractical for a large permutation set, due to excessive storage demand. To avoid this, we propose the usage of an one-to-one mapping method, which is implemented by a Heap’s-like Algorithm.
