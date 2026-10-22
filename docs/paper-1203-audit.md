# Audit of Schuch et al., arXiv:1203.4816v3

Source: the uploaded **Resonating valence bond states in the PEPS formalism**, Norbert Schuch, Didier Poilblanc, J. Ignacio Cirac, and David Pérez-García, 20 PDF pages. Page numbers below are PDF page numbers, which agree with the printed page numbers. This document audits the implementation present when the papers were supplied; subsequent new calculations must be documented separately.

## What the paper establishes and what it actually plots

The paper focuses on **kagome**, not a comparative numerical study of square and kagome correlation exponents. The square lattice also appears as the computational lattice obtained by blocking kagome, and as a lattice for the toric code. Those mentions must not be interpreted as square-lattice RVB correlation calculations.

On PDF p.7, Sec. IV B 1 says: “We have computed various correlation functions, all of which show an exponential decay with a slope which changes smoothly with θ”. This supports qualitative exponential decay throughout the kagome interpolation. It does not supply a universal numerical dimer correlation length or the values previously fitted by our program.

**Figure 11, PDF p.7:** the text identifies “the connected dimer-dimer correlation function between two z ⊗ z operators [z = diag(1,−1,1)]”, with each dimer supported on one blocked square tensor of Fig.10b and the two squares in the same row. The plot is for the infinite system, contracted using an infinite MPS and iTEBD, with boundary bond dimensions **D=40 and D=60**. Curves run from θ=0 to θ=0.9. The x axis is distance **between blocked tensors**, not simply distance between individual kagome spins. This contraction D is not the RVB PEPS virtual dimension, which is 3.

The original adapter `scripts/run_upstream.py` measures the physical scalar bond operator **S_i·S_j**. These curves therefore **are not a direct reproduction of Fig.11**. They can independently support the qualitative short-range conclusion, provided geometry, units, fit interval and approximation are reported.

The paper's three-level z and Eq.(12)'s four-dimensional/two-qubit local physical space require an explicit embedding before claiming an exact operator match. In the original qutrit-pair basis, `z ⊗ z` on one vertex acts as physical Pauli Z on both |s,2> and |2,s>. In contrast, z on just one qutrit is a location-dependent operator. The sentence describing Fig.11 does not specify a bond label or the full embedding. A new physical Z_i Z_j calculation should state its convention and not silently identify it with a uniquely reconstructed paper operator.

Even at the singlet endpoint, a four-spin ZZ connected function cannot in general be obtained by multiplying the scalar S·S connected function by 16/9. The two-spin relation <Z_i Z_j>=(4/3)<S_i·S_j> follows from SU(2), but the four-spin function also contains the rank-two component of S_i^z S_j^z. Its contribution can have a different decay. It requires an explicit insertion and contraction.

The PDF contains the plotted curves as vector paths. `scripts/extract_paper_fig11.py` reproducibly extracts their displayed vertices into `results/paper-recheck/published_fig11.json`, including the PDF checksum, raw coordinates and axis calibration. These are **published plotted coordinates, not the authors' original raw data**. At θ=0 the D=60 curve gives approximately C(1)=0.0658, C(2)=8.84×10^-4, C(3)=2.43×10^-5, C(4)=1.52×10^-6. D=40 gives approximately C(4)=1.87×10^-6, so the longest displayed distance is itself sensitive to the contraction cutoff. Points clipped below 10^-6 are absent from the extracted dataset, not reported as zero. The script validates the ten curves of each color and their theta order from steepest (0.9) to shallowest (0).

**Figure 12, PDF p.8:** the caption calls the plotted quantity a transfer gap Δ, but explicitly defines it as “the ratio of the third eigenvalue and the average of the two maximal eigenvalues.” The two leading eigenvalues arise from parity/topological sectors and are nearly degenerate for even circumference. The curves use N_h=4,…,10; the infinite-circumference fit uses N_h=6,…,10 and a+b exp(−c N_h). Thus the figure must not be compared with an unqualified λ_2/λ_1 from another transfer geometry. To turn a decaying eigenvalue ratio ρ into a length use ξ=−a/log|ρ|, with step a and the observable's symmetry channel specified. Neither the plot's Δ nor 1−Δ is automatically a measured dimer length.

## State and tensor conventions

On PDF pp.2–3, Eqs.(1) and (2), the elementary triangle state and vertex projector are

\[
|\epsilon\rangle=\sum_{ijk}\epsilon_{ijk}|ijk\rangle+|222\rangle,
\qquad P=|0\rangle(\langle02|+\langle20|)+|1\rangle(\langle12|+\langle21|).
\]

Here ε_012=1; |2> means no singlet on that virtual leg, **not spin S_z=0**. These are the tensors constructed in `rvb/peps.py` by `triangle_tensor()` and `kagome_site_tensor()`. The implementation blocks three physical spins and two triangle tensors, as the paper does in Fig.10, although axes, selected bond and distance direction still need matching before pointwise comparison.

The paper orders each triangle clockwise (p.3, discussion of Eq.(1)); `rvb/lattices.py` orders every triangle counterclockwise. This alone does **not** change the normalized RVB state: all singlet orientations are reversed, and every N-spin covering has N/2 singlets, so every amplitude receives the same overall factor (−1)^(N/2). Flipping only a subset of arrows would require a separate gauge/sign check.

On PDF p.3, Eq.(4) and Sec. II D explain that normalizing all singlets changes the total state only by an overall constant because every covering contains the same number of singlets. The exact covering oracle's normalized singlets and the PEPS's unnormalized singlets are therefore compatible for normalized observables.

## Orthogonal dimers and the interpolation

On PDF p.2, Sec. II A defines the dimer state using orthogonal covering states |D>, and defines the RVB state as the equal-amplitude sum of nonorthogonal singlet products |σ(D)>.

On PDF p.3, Eq.(3), P_perp keeps the local states |02>, |12>, |20>, |21>, rather than erasing the triangle/location information. The text explicitly stresses that this enlarged-space dimer state still “uses actual singlets built of two entangled particles.” Thus orthogonalizing coverings does **not** mean replacing a spin singlet by a product state or assuming the original spin-1/2 singlets suddenly become orthogonal.

For a diagonal occupation operator n_e, orthogonality eliminates cross-covering interference, and the equal-amplitude state has the uniform classical dimer occupation distribution. This justifies `rvb/orthogonal.py` as an occupation benchmark. It does not establish equality of occupation and scalar-spin correlation amplitudes between different Hilbert-space representations. Distant zero occupation correlations on kagome do not imply absence of entanglement or of all local correlations.

On PDF p.6, Eq.(12), the interpolation is

\[
P(\theta)=|+\rangle\sum_{s=0,1}|s\rangle(\langle s2|+\langle2s|)
+\theta|-\rangle\sum_{s=0,1}|s\rangle(\langle s2|-\langle2s|).
\]

**θ=0 is RVB; θ=1 is the orthogonal dimer endpoint**, up to local isometries. An implementation must retain the additional physical label for intermediate θ. Varying the amplitude of ε_ijk versus |222> is not this interpolation. Nor does running a separate classical endpoint solve the intermediate-θ PEPS.

There is an apparent limitation or mismatch in the prose overlap shorthand immediately following Eq.(12), on PDF p.7. Read literally, it proposes covering overlap x^(L−2N) with x=(1/√2)(1−θ²)/(1+θ²). Equation (12) itself assigns normalized flag states (|+>±θ|−>)/√(1+θ²) to the two possible occupied triangles at each vertex. Their overlap is q=(1−θ²)/(1+θ²). Thus two fixed coverings acquire q^m, where m counts vertices at which the occupied triangle changes, multiplying their original spin-singlet overlap. This is an exact local consequence of Eq.(12).

For an alternating six-edge loop whose occupied triangle changes at all six vertices, the normalized spin overlap is 1/4 and the Eq.(12) result is q^6/4. The literal scalar-x shorthand instead gives q^4/4. At θ=0.5 these are 0.011664 and 0.0324. An independent dense six-spin plus six-flag construction gives 0.011664. Such six-flag covering differences also occur on an actual 3×2 kagome torus; exact covering sums give the same spin overlap 1/4. Conversely, the six-edge loops in the smallest 2×2 kagome torus can change only four flags, so testing only that cluster can mask the distinction. This audit therefore follows **the explicit tensor Eq.(12)** and does not force its output to agree with the prose shorthand. Establishing whether the prose intended a more restrictive loop class or another parametrization would require further source/author clarification.

## Units and boundaries

The implementation uses kagome Bravais vectors a_1=(2,0) and a_2=(1,√3) when nearest-neighbor spacing is 1. Each elementary axial unit-cell translation therefore has length **two nearest-neighbor spacings**. Its quoted ξ_x≈0.624 cells corresponds to ≈1.248 nearest-neighbor spacings along that chosen axis, but remains a finite-window estimate for the program's chosen scalar bond, not a paper value.

The paper's cylinder circumference N_h counts blocked sites in Fig.10. It uses left/right periodic connections and upper/lower virtual |2> boundary conditions before taking the infinite-length limit (p.7). Our exact cylinder calculations have finite circumference and an infinite axial direction, whereas the existing CTMRG plots approximate the infinite plane. Neither a circumference label nor a contraction bond dimension is interchangeable with the other.

The finite-width exact-cylinder calculation is useful validation, but its exponential tail cannot by itself establish an infinite-plane asymptotic law. Likewise, the finite-χ square CTMRG cutoff is not evidence of a finite physical bulk dimer length.

## Changes required before a quantitative reproduction claim

1. Keep the original scalar S·S and classical occupation results explicitly labeled as separate observables and models.
2. Contract the four-spin observable explicitly, with local embedding, blocked bond and translation direction recorded.
3. Implement Eq.(12) if comparing the θ-dependent curves; verify both endpoints and the state convention independently.
4. Match Fig.10 geometry and compare at least two contraction dimensions. A plot of matching qualitative decay is insufficient to establish pointwise reproduction.
5. If comparing Fig.12, treat the two leading topological/parity eigenvalues according to its caption and match the number of transfer steps and circumference before extracting lengths.
6. Replace the old statement that primary paper access remains blocked: these two user-supplied PDFs are available and have now been read. Additional papers may still be needed for other numerical benchmarks.

The present paper is sufficient to specify the kagome tensors, endpoint distinction and interpolation. For the exact short-range classical kagome result, its reference 11 is Misguich, Serban and Pasquier, Phys. Rev. Lett. **89**, 137202 (2002), arXiv **cond-mat/0204428**; that is an optional useful original-source check, not a blocker for independently verified occupation calculations. Square-lattice exponents need the square-specific sources cited in the companion paper rather than attribution to this kagome-focused work.
