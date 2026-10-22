# Exact cylinder checks of the square-octagon singlet RVB

These calculations use the requested conventional equal-coefficient oriented-singlet state. The raw physical tensor is contracted exactly. There is no boundary-bond truncation, no semionic loop sign, and no inserted vison. All lengths below are finite-cylinder channel lengths, in original square-lattice unit cells. They are not a bulk extrapolation.

The cylinder wraps the y direction. Even and odd refer to virtual-vacuum parity, separately fixed to the same value in ket and bra. Distances run along x. Local bond operators are scalar S_i dot S_j, with ring-ring sites (0,3), and ring-dummy sites (0,4) or (3,4). Spin lengths use S_z at ring site0.

| Width | Parity | Ring-ring mean | Bond ξ | Ring-spin ξ | Dominant residual |
|---:|---|---:|---:|---:|---:|
| 2 | even | -0.329918524 | 0.906704390 | 0.217571563 | 1.1e-16 |
| 2 | odd | -0.250000000 | 0.835945134 | 0.224499783 | 2.5e-16 |
| 3 | even | -0.319061568 | 0.789440305 | 0.210311074 | 4.1e-16 |
| 3 | odd | -0.268574419 | 0.593664999 | 0.213154315 | 1.2e-15 |
| 4 | even | -0.307434131 | 0.834939009 | 0.229809448 | 8.5e-16 |
| 4 | odd | -0.284795250 | 0.622156688 | 0.231518112 | 2e-15 |

All three tested bond operators share the leading coupled scalar eigenmode at each width and sector. Spin correlations alternate in sign and have much shorter lengths. The full sector eigenspectrum is used, including both left and right modes and degenerate eigenspaces. A generic second transfer eigenvalue would mix physical channels and can be dominated by sector splitting.

For the width2 cylinder, the selected dummy-site S_z autocorrelation has no resolved nonzero channel; the output records ξ as null rather than inferring a length from roundoff. At width4 odd parity the dummy-site spin couples to a different leading triplet mode than the ring-site spin, giving ξ=0.223780862 rather than0.231518112.

Eigenmode reconstructions agree with directly propagated centered correlations to better than3e-15 for separations r>=2. Dominant fixed-point residuals are below2.1e-15. The exported JSON contains r=1..32, individual coupled eigenvalue ratios, spectral weights, and left/right mode residuals. Plot only values above the stated1e-13 absolute floor; tiny values below it are not evidence for a long tail.

The even/odd and circumference dependence remains material. The results support short-ranged physical channels on the cylinders tested, but do not establish a unique thermodynamic correlation length.

Validation:9 tests in tests/test_octagon_transfer.py pass, including independent covering-sum norms and bond means in each sector on1x1 and1x2 tori, alongside eigenmode reconstruction of physical correlations.
