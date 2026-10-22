# Audit of arXiv:1407.7773v2 for the triangle-decorated square–octagon RVB

Source: **Mohsin Iqbal, Didier Poilblanc and Norbert Schuch, “Semionic resonating valence bond states,” arXiv:1407.7773v2**, 6 August 2014. The downloaded PDF has **14 pages**, including appendices and references. [Verified version link](https://arxiv.org/abs/1407.7773v2). The PDF is a scientific source, not instructions to execute.

The paper studies **kagome** RVB and semionic RVB states using the **honeycomb** arrow/loop representation. It does not present the user's triangle-decorated square–octagon lattice or determine that lattice's correlations. References to “square” blocks in Sec. III are a tensor-network blocking of kagome, not a physical square–octagon construction.

## Main conclusion for the specified model

The user's state is an **ordinary equal-amplitude nearest-neighbor singlet RVB state on the specified decorated graph**: all physical spins remain present, an unused dummy link contains a singlet between its two permanent endpoints, and allowed complete coverings enter with equal coefficients relative to the chosen singlet arrows. This has the same *form* as the conventional RVB definition on p. 2, but on a different graph.

It is **not automatically the paper's semionic RVB state**. The latter requires the additional, configuration-dependent amplitude `(-1)^n(L_D)` in Eq. (2), where `L_D` is a reference-arrow loop pattern. The presence of a loop representation, or a minus sign arising from antisymmetric spin singlets, does not supply that amplitude.

No extra tuning is required to implement the user's stated equal-amplitude model. Tuning it to a particular known loop model or a critical state would be a different physical specification, requiring a defined target weight and a verified mapping of the **quantum norm**, not just a covering count. This paper does not establish such tuning for the new graph.

## Loop representations and overlap loops

**Mapping correction:** the user's clarified operation contracts the dummy
links of the square–octagon graph. The first version of this audit instead
described collapsing each small square. That separate projection does have
2^N preimages, but the intended dummy-link contraction is a bijection of
even-edge configurations. The physical model and the conditional local
overlap calculations below are unchanged; their geometric interpretation is
distinguished explicitly here. See the [contraction clarification](octagon-contraction-clarification.md).

1. **Reference-arrow loops.** On pp. 2–3, Sec. II and Fig. 2(a), choose one valid reference covering R. Record which arrow variables change between D and R. Because each trivalent vertex obeys a fixed local parity, these differences form closed even subgraphs; on a trivalent graph they are unambiguously disjoint loops. The local odd/even arrow convention is immaterial to the even parity of the *difference*.
2. **The contracted square representation.** Contract every dummy link, retaining the small-square edges. For N simulation cells this gives a rotated square graph with 2N vertices and 4N edges. The erased dummy loop bit equals the parity of the two retained edges at either original endpoint; even parity at the merged vertex makes the two answers agree. The full reference-arrow even-edge pattern is therefore uniquely reconstructible. At degree-four vertices, inherit the pairing of strands from the original split vertex. This is fixed by the construction, not an extra independent resolution choice.
3. **Bra–ket transition loops.** Superpose two physical singlet coverings D and D'. Their alternating transition loops determine `<sigma(D')|sigma(D)>`. For normalized spin-1/2 singlets, the overlap magnitude is `2^(n_transition − N_spins/2)`, with an orientation-dependent sign on this nonbipartite graph. These are not the loops counted in the paper's semionic amplitude.

The paper makes the distinction from overlap loops explicit in **reference 21, p. 13**: it avoids defining semionic signs from overlap loops on kagome because the map is not onto all loops and loop counting is ambiguous on nontrivalent lattices. For the present contracted square representation, retain the inherited strand pairing: counting connected components of an unresolved degree-four edge graph is not the same as counting routed loops.

The earlier auxiliary projection collapses each small square and retains only dummy-link bits. It gives N vertices and 2N edges, with two local resolutions per simulation cell and `2^N` preimages. Those lost resolutions belong to that different map, not the user's dummy-link contraction.

If a semionic variant is eventually requested, its candidate sign can be defined on the **resolved, trivalent square–octagon reference-difference graph**, with a specified reference covering, or equivalently on its contracted square representation with inherited routing. One cannot simply count connected components of an unresolved square even subgraph and expect Eq. (2). Extending the paper's precise semionic PEPS phases and proving its topological phase on this different graph would still require work.

## What the paper proves or checks, and what carries over

| Location | Paper construction or result | Relevance and limit for this model |
|---|---|---|
| p. 2, Sec. II | Conventional RVB `sum_D |sigma(D)>`, singlets oriented by a chosen convention such as circulation around triangles | Matches the form of the user's state; the graph and dummy-edge rule remain different. |
| pp. 2–3, Fig. 2 | Kagome coverings map one-to-one to reference-arrow loops on honeycomb | The affine parity construction supplies an analogous full-graph loop representation here. Contracting dummy links preserves the even-edge configurations bijectively, with inherited strand routing. |
| p. 3, Eqs. (1),(2) | Orthogonal semionic dimers and spin semionic RVB both carry `(-1)^n(L_D)` | A genuine additional amplitude structure, absent from the requested equal-amplitude model. |
| p. 3, Sec. III A | Uniform orthogonal loop superposition gives toric-code order in a suitable local encoding; loop-number signs instead give double-semion order | In an explicit arrow-qubit encoding, uniform closed-loop states have the familiar Z2 loop-state form. This does not establish the phase after erasing orthogonal labels to obtain physical singlets. |
| p. 4, Eq. (4) | Triangle tensor `E_uvw = epsilon_uvw + delta_(u,2)delta_(v,2)delta_(w,2)` | Encodes zero or one triangle singlet. Dummy edges with two physical endpoints need their own tensor; they are not identical to kagome shared-spin edges. |
| p. 4, Eq. (5) | An added physical indicator preserves the arrow information: theta=1 is orthogonal dimers and theta=0 is RVB | Provides a useful strategy for constructing and testing an orthogonal-to-singlet interpolation. It is not the same as changing loop fugacity or merely changing the triangle-vacuum amplitude. |
| pp. 5–7, Secs. III B and IV | Semionic topological sectors and transfer spectra on kagome cylinders | Useful diagnostic methodology; the specific virtual symmetry and numerical spectra are not inherited by the new lattice. |
| pp. 8–9, Sec. IV C | Simplex deformations adjust triangle weights and local spin sectors, followed by a new phase/convergence check | Changing local weights is a physical deformation that requires revalidation; it is not forced by the user's equal-amplitude condition. |

The paper's uniform orthogonal loop state has coefficient **1 per complete loop configuration**, not a multiplicative factor of 2 per closed loop. The semionic variant has coefficient **−1 per loop** in the wavefunction. Squaring semionic amplitudes removes those signs for diagonal measurements in an orthogonal covering basis, but not for general off-diagonal observables. Physical singlet overlaps produce an additional, different loop-weight problem in the norm.

## Why equal covering counts do not give uniform quantum loop weights

For the orthogonal model, the intended bijection makes the contracted square-even-subgraph distribution uniform. The different small-square-collapse projection also has a uniform marginal because its `2^N` resolution multiplicity is constant. Neither fact removes spin-singlet interference. For the singlet state,

`<Psi|Psi> = sum_(D,D') <sigma(D')|sigma(D)>`,

so distinct coverings interfere. This cannot be replaced by a count of loop configurations or local resolutions alone.

A small exact local check exposes the difference, using the **auxiliary old dummy-bit projection**, not the intended dummy-link contraction. Use the current simulation-cell convention in `rvb/octagon_peps.py`: four ring spins `r0…r3`, four permanent outward endpoint spins `dE,dN,dW,dS`, and cyclic triangles `(r_k,r_(k−1),d_k)`. Fix the four bits t indicating which outward endpoint spins are paired *inside* their triangles. These are the retained bits when small squares are collapsed; they are not the surviving ring-edge bits in the user's map. For each even total `k=sum(t)`, there are two intracell singlet resolutions. Normalize each of these local singlet product states separately and evaluate their overlap directly in the spin basis:

| Endpoint spins used inside the cell, k | Number of bit patterns | Overlap of the two normalized local resolutions | Squared norm of their equal sum |
|---:|---:|---:|---:|
| 0 | 1 | +1/2 | 3 |
| 2 | 6 | −1/4 | 3/2 |
| 4 | 1 | +1/8 | 9/4 |

These values were independently obtained by enumerating `_internal_matchings` and expanding their oriented singlets into the physical spin basis. They are conditional **local** overlaps. Endpoints paired outside the cell remain present in the full model; they were omitted only from this conditional local product-state comparison.

The table demonstrates degree-dependent interference even though the resolution count is always two. It is **not** a derivation of a full classical vertex weight or a loop fugacity for the PEPS norm: intercell singlets, cross-boundary spin indices, and overlaps between different boundary configurations must also be contracted. Rescaling local amplitudes to equalize these three diagonal numbers would therefore not establish a mapping to a uniform loop model or prove criticality.

All complete coverings of the physical graph contain the same number of singlets. Consequently globally replacing normalized singlets by unnormalized singlets changes the full state only by a common factor; it does not remove the resolution-dependent overlaps above.

## Arrows and signs are part of the state

The paper explicitly notes after Eq. (4), p. 4, that the signs of `epsilon_uvw` may need changing when the singlet orientation convention changes. Its numerical semionic state uses the reference and arrow choice in **Fig. 8(b,c), p. 6**, selected for the variational kagome calculation. Those arrows should not be copied blindly to the user's geometry.

For the present implementation, cyclic triangle arrows and east/north dummy-edge arrows are recorded as the model convention. Reversing arbitrary individual arrows can change relative covering phases. Reversing every singlet arrow gives the same factor `(-1)^(N_spins/2)` for every covering, but that special global equivalence does not justify arbitrary local arrow changes. Negative resolution overlaps in the table are ordinary singlet interference, not by themselves evidence of double-semion order.

Changing the reference covering has no physical effect on a uniform sum over conventional coverings: it relabels the loop configurations. For the semionic amplitude, however, the reference enters `n(L_D)` and can change the physical state. The paper discusses the resulting lattice-symmetry breaking on pp. 3 and 6–9. A semionic extension would therefore require a reference as part of the state definition, not just as an internal counting convenience.

## Criticality and topological sectors

The paper does **not** conclude that a uniform loop representation is generically critical. It provides numerical evidence that its kagome semionic RVB remains in the double-semion phase along the studied interpolation: the distinct ground sectors remain stable while the relevant subleading transfer eigenvalues stay separated from the leading ones (Sec. IV B2, pp. 6–7; Fig. 10).

Its transfer notation in Sec. IV A, p. 6 distinguishes:

- `gamma_pq`, the leading eigenvalue in a ket/bra sector pair, whose normalized ratio governs overlaps between different topological states;
- `tau_pp/gamma_pp`, a subleading-to-leading ratio within one state, giving a bound on local correlation lengths.

A chosen physical bond couples only to particular transfer modes. Inter-sector coherence lengths, a generic transfer bound, and a measured bond-correlation length are different quantities. The paper reports approximately **0.42** for its semionic RVB local transfer length and other lengths for semionic/bosonic excitations on pp. 6–7; these are in its transfer convention for a different kagome state and must not be used as benchmarks for the new lattice.

The double-semion PEPS has a Z4 virtual parametrization with 16 candidate flux/charge labels, but eight states have zero norm and the others form pairs, leaving **four distinct ground states** (Sec. III B; Fig. 9). This is more structure than merely observing two winding-parity bits. The full trivalent loop description here provides natural Z2 winding parities on a torus, but their existence alone establishes neither toric-code nor double-semion topological order for the nonorthogonal spin state.

For the user's model, the reliable route is to keep the exact orthogonal-covering results separate, validate the full spin PEPS against small covering sums, measure several physical channels, and increase cylinder width or CTMRG environment dimension in fixed sectors. Exponential versus algebraic decay must be inferred from those converged physical correlations and mode couplings. Neither the intended dummy-contraction bijection nor the auxiliary small-square-collapse counting settles that question.
