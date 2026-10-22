# Exploratory three-dimensional spin-one overlap sampling

These are actual calculations on periodic cubic boxes, with 3 L^3 spin-one
sites. They are not classical covering correlations and do not use a narrow
transfer cylinder. They also do **not** establish a thermodynamic correlation
length: long-distance signals on L=6 are statistically unresolved, and rare
replica sectors have not passed a full mixing study.

The state is the equal-coefficient sum of the allowed normalized coverings,
with Cartesian dimer delta/sqrt(3), trimer epsilon/sqrt(6), and fixed x,y,z
trimer ordering. Each individual covering-pair contraction is exact integer
arithmetic before restoring the square-root normalization.

`summary.json` contains spin correlations and 20-block delete-one-block
jackknife errors after discarding the first quarter of each run. These are
statistical error estimates, not bounds on equilibration bias. Each measurement
uses a random origin and edge axis; the displacement is parallel or perpendicular
to that axis. Cubic symmetry makes this equivalent to averaging the x-edge/x-axis
correlator over equivalent directions. Spin units are hbar=1 and the measured
quantity is <S_i dot S_j>.

For biased runs, the sampling target is |<bra|ket>| B, with B=100 if the two
coverings have different trimer-location fields, and B=1 otherwise. The recorded
spin estimator is <bra|S_i dot S_j|ket>/|<bra|ket>|. The physical ratio is
sum(weight * estimator) / sum(weight * sign), where weight=1/B. Trimer densities
and replica differences are covering-space mixing diagnostics, not physical
quantum observables.

Reproduce a principal run from the repository root:

```sh
python scripts/octahedral_bulk_overlap_sampling.py --L 4 --steps 100000 --initial diagonal --seed 39383
python scripts/octahedral_bulk_overlap_sampling.py --L 4 --steps 100000 --initial independent --seed 49383
python scripts/octahedral_bulk_overlap_sampling.py --L 6 --steps 20000 --initial diagonal --seed 49313
python scripts/octahedral_bulk_overlap_sampling.py --L 6 --steps 30000 --initial diagonal --seed 99295 --bias 100
```

The original files without a `_bias1` suffix predate that naming convention.
The physical estimator was subsequently rewritten in logarithms to avoid
underflow; this does not change these finite-box results.

## Sign and contraction checks

- 430 independently sampled allowed covering pairs on L=3 and 30 on L=4 had
  positive exact overlaps. No sign-free theorem is claimed.
- For L=3, 400-pair mean contraction time was 0.090 s; maximum intermediate
  had 2187 entries. For L=4, 30-pair mean was 0.356 s and maximum intermediate
  had 59049 entries. Independent overlaps are exponentially tiny and are
  unsuitable for a naive independent-proposal estimator.
- A coordinate-line flip of a diagonal covering has exact overlap
  3^[1-(L+t)/2], where t is its number of trimer centers. These winding
  contributions can cause substantial finite-size effects.
- A cube flip has exact overlap (1/54) 3^(-n/2), where n counts incoming
  external ports at the cube corners. Its mean over uniformly sampled
  diagonal backgrounds is 6.212007830959485e-5 for L>=3.
- Joint and single replica line/cube moves span the full covering space for
  L=3 and L=4, checked by exact binary rank. Irreducibility does not establish
  practical mixing.
- Unit tests check diagonal norms, dimer and trimer spin correlations, raw
  unreduced tensor contractions with spin insertions, Hermiticity, the known
  cube overlap, and the line overlap formula.
