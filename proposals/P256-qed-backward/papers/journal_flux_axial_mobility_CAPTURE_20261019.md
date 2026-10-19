# Intermediate flux–axial–mobility paper capture

This is a byte-preserved public capture of **intermediate finite-model research** associated with [issue #222](https://github.com/vantasnerdan/substrate-framework/issues/222). Issue #222 remains **ACTIVE**: its actual coupled physical goal—an explanatory, quantitatively related interacting charged spin-half matter/photon sector—is unmet. This capture is not completion of that goal, a derivation of nature or continuum QCD, or a journal submission or acceptance.

## Read and reproduce

- [Manuscript PDF](journal_flux_axial_mobility_20261019/journal_flux_axial_mobility_20261019.pdf)
- [Manuscript TeX](journal_flux_axial_mobility_20261019/journal_flux_axial_mobility_20261019.tex)
- [Self-contained sealed ZIP](journal_flux_axial_mobility_20261019.zip)
- [Reproduction instructions and environment requirements](journal_flux_axial_mobility_20261019/journal_flux_axial_mobility_REPRODUCE_20261019.md)
- [Existing fresh-environment reproduction evidence](journal_flux_axial_mobility_reproduction_20261019.json)
- [Payload manifest](journal_flux_axial_mobility_20261019/journal_flux_axial_mobility_manifest_20261019.json)

The ZIP SHA256 is `331ca5663a07af6558d64df153a211764a2602e9ef41bd2e16c8dbc1bbf4fe0b`. Its exact 21 members (20 payload files plus manifest) are also extracted alongside it. No recursive research-directory copy or scientific-source rewriting was used. The reproduction evidence is a separate byte-preserved sidecar, not an additional archive member.

From this directory, verify the frozen payload and selected scientific source/result/figure bindings without repeating the calculations:

```sh
python journal_flux_axial_mobility_20261019/journal_flux_axial_mobility_package_20261019.py verify
```

The capture-time invocation reported `20` verified payload files and matched source bindings. The separate existing reproduction JSON records an actual fresh-environment `all` execution: three scientific providers, figure/table generation, and TeX/BibTeX compilation; scientific comparisons passed at relative tolerance `2e-7` and absolute tolerance `2e-8`, with frozen references unchanged. Publication did not rerun that full calculation. Follow the linked instructions to perform a new reproduction; regenerated PDF bytes are not promised bit-identical.

## What the result does and does not establish

The paper specifies a canonical SU(3) bouquet–tree Hamiltonian with twelve positive-energy CAR modes per matter root and an adopted nonrelativistic Pauli vertex. Its contribution connects explicit chromomagnetic insertions in physical transfer resolvents to baryon mobility, axial deficit and local splitting, then calculates canonical higher-flux local axial response and a compact six-quark versus separated proton–neutron comparison using the same local coupling. Baryon hopping itself is established prior work; no priority claim for that mechanism is made here.

The transport expansion and nonlinear local observables are distinct consumers of the specified law, not a solved finite-hopping joined spectrum. Canonical two-source results are finite Galerkin calculations; the manuscript distinguishes them from claims about an untruncated canonical tower. An externally supplied constituent mass and effective Pauli interaction are not derived Dirac-sea matter dynamics. The charged matter/photon feedback, continuum matching and broader physical explanation required by issue #222 remain open.

Only the manifest-selected scientific payload, sealed archive and reproduction evidence are captured here. Opaque provenance identifiers within scientific sources/results are preserved as nonsecret identifiers. Private notes, memory/worldmodel material, campaign narratives, credentials and unrelated research outputs are not part of this publication.
