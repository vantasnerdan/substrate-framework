# Reproducing the flux–axial–mobility paper

This package contains the anonymous manuscript, final PDF, bibliography, numerical figure/table generator, frozen calculation results and the scientific Python sources needed to reproduce them. It is a local, explicit-file archive, not an installer. Extracting it does not execute code, install packages, download references or contact services. No public upload is performed by the entrypoint.

## Supported reference environment

The observed scientific environment is **CPython 3.12.3**, NumPy **2.5.3**, SciPy **1.18.1**, threadpoolctl **3.7.0** and Matplotlib **3.11.2**, on Linux x86-64. These are measured installed versions, not inferred minimum requirements. The supplied requirements file pins these four direct dependencies; it is not a lockfile for their transitive dependencies, native BLAS libraries or operating system. Other Python/library combinations are unverified, not claimed unsupported in principle.

The observed paper toolchain is **pdfTeX 3.141592653-2.6-1.40.25** and **BibTeX 0.99d**, TeX Live 2023/Debian. Rebuilding requires `pdflatex` and `bibtex` on PATH, the standard `article` class, `geometry`, `amsmath`, `amssymb`, `bm`, `booktabs`, `graphicx`, `natbib`, `hyperref`, and the `unsrtnat` bibliography style. REVTeX and a network bibliography lookup are not required. A complete TeX distribution is a separate prerequisite; it is not vendored in the ZIP.

For an explicit Python environment installation, from a directory outside the extracted bundle:

```sh
python3.12 -m venv journal-venv
. journal-venv/bin/activate
python -m pip install -r extracted-paper/journal_flux_axial_mobility_20261019/journal_flux_axial_mobility_requirements_20261019.txt
```

Environment creation and `pip install` are deliberate user commands, not unpack hooks. The builder and hash verifier themselves use only the Python standard library. The reproduction command checks installed distribution metadata and tool availability before starting; the subprocesses use the same Python interpreter as the entrypoint (`sys.executable`), with argv lists and no shell. BLAS/OpenMP child thread counts are set to one. Matplotlib uses its noninteractive Agg backend.

## Assemble a complete archive (maintainer)

Run these commands from the directory containing the paper and package sources:

```sh
python journal_flux_axial_mobility_package_20261019.py status
python journal_flux_axial_mobility_package_20261019.py build
```

`status` reports file-presence prerequisites and writes nothing; it exits nonzero when any required final file is missing. `build` additionally checks scientific source/result and figure-data provenance bindings. In particular, it refuses to produce a partial archive when the final manuscript PDF is absent. If the bibliography, scientific source, result JSON or figure generator changes, regenerate the dependent figure/table/data outputs before building. The final PDF must be assembled from those final inputs by the manuscript integration owner; the builder does not claim to establish its compilation history merely by hashing it.

A successful build writes exactly:

- `journal_flux_axial_mobility_20261019.zip`, containing one portable directory named `journal_flux_axial_mobility_20261019`;
- `journal_flux_axial_mobility_manifest_20261019.json`, also included inside that directory.

The build prints the ZIP SHA256 and member count. `--output another-name.zip` changes only the ZIP destination; the manifest sidecar remains beside the entrypoint. Inputs are read from beside the entrypoint, independently of the invocation directory. Existing designated ZIP/manifest outputs are replaced only after all input prerequisites and bindings pass. No recursive research-directory copy is used.

## Extract and verify frozen references

From the directory holding the ZIP:

```sh
python -m zipfile -e journal_flux_axial_mobility_20261019.zip extracted-paper
python extracted-paper/journal_flux_axial_mobility_20261019/journal_flux_axial_mobility_package_20261019.py verify
```

The JSON manifest records SHA256 and byte count for **every payload file**, including the entrypoint, these instructions, requirements, scientific source/results, bibliography, manuscript and PDF. The manifest excludes itself and the enclosing ZIP to avoid circular hashes. The printed outer ZIP digest can be retained separately. These are integrity checks, not a signature or proof of scientific correctness. `verify` compares the frozen bytes and the scientific source/result/figure provenance; it does not rerun calculations.

## Full calculation-to-paper reproduction

With the reference environment activated and TeX installed, run once from outside the frozen bundle:

```sh
python extracted-paper/journal_flux_axial_mobility_20261019/journal_flux_axial_mobility_package_20261019.py reproduce --mode all --work-dir paper-rerun
```

`paper-rerun` must not already exist, and must be separate from (not inside or an ancestor of) the frozen directory. This protects both the archived reference and earlier reruns. The command:

1. Verifies the complete frozen manifest and provenance bindings.
2. Copies only the explicitly named scientific and manuscript inputs into the new work directory.
3. Executes `rich_journal_flux_carrier_20261019.py`, then `rich_journal_baryon_band_20261019.py`, then `rich_journal_compact_hierarchy_20261019.py`, producing three fresh calculation JSONs. The compact calculation reads the newly produced carrier JSON for its mass-three nucleon baseline.
4. Compares all 24 carrier rows' energy, source energy, elastic axial vertex, vacuum weight, probability-law defect, complete available inelastic axial response and degree-three exterior diagnostics, plus the complete band-result object, against the frozen references. For the compact hierarchy it compares both masses' compact six-quark energy, single three-quark cell energy, separated proton–neutron energy, compact-minus-separated penalty, one-source and leading-order penalties, and physical-sign label. Numeric comparison uses relative tolerance `2e-7` and absolute tolerance `2e-8`; row identity and nonnumeric labels must match exactly. Solver timings, matrix-vector call counts, carrier runtime metadata and compact provenance metadata are not scientific comparison targets. In particular, fresh carrier-result timings legitimately change the compact result's input hash without changing merger physics. Any scientific mismatch is recorded and stops the run before figure generation.
5. Executes `journal_flux_axial_mobility_figures_20261019.py` on the fresh results, recreating the figures PDF, numerical data JSON and TeX tables.
6. Executes `pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error`, then `bibtex`, then two further identical `pdflatex` passes on `journal_flux_axial_mobility_20261019.tex`.
7. Rechecks the frozen reference and writes `paper-rerun/reproduction_report.json` with the actual environment, scientific comparison outcome, output SHA256s and success/failure status.

Subprocess stdout and stderr are captured separately as `*.stdout.txt` and `*.stderr.txt` in `paper-rerun`; TeX also writes its normal auxiliary files there. The final recreated paper is `paper-rerun/journal_flux_axial_mobility_20261019.pdf`. Frozen references are never overwritten. Fresh source-bound data hashes legitimately differ when regenerated scientific result bytes contain new timings. Figure/PDF bytes can differ because of timestamps, fonts and TeX/Matplotlib versions; bitwise reproduction of generated PDFs is not asserted. The scientific programs perform actual finite-Hamiltonian calculations rather than returning cached reference data.

For figures/tables/paper only, using the frozen scientific results without repeating the calculations:

```sh
python extracted-paper/journal_flux_axial_mobility_20261019/journal_flux_axial_mobility_package_20261019.py reproduce --mode paper --work-dir paper-from-reference
```

This mode explicitly reports that calculations and scientific numeric comparison were not executed. It still regenerates figures, tables, data, bibliography and manuscript PDF in a new directory. A successful PDF process is not a substitute for reviewing the mathematical argument and scientific limitations in the manuscript.

## Exact archive contents and dependency boundary

The archive builder uses this finite allowlist, plus the generated manifest:

```text
journal_flux_axial_mobility_package_20261019.py
journal_flux_axial_mobility_REPRODUCE_20261019.md
journal_flux_axial_mobility_requirements_20261019.txt
rich_journal_flux_carrier_20261019.py
rich_journal_baryon_band_20261019.py
rich_journal_compact_hierarchy_20261019.py
rich_nucleon_hierarchy_20261019.py
rich_composite_local_flux_binding_20261017_model.py
rich_colour_spin_root_20261019.py
rich_journal_flux_carrier_20261019_result.json
rich_journal_baryon_band_20261019_result.json
rich_journal_compact_hierarchy_20261019_result.json
journal_flux_axial_mobility_figures_20261019.py
journal_flux_axial_mobility_figures_20261019.pdf
journal_flux_axial_mobility_data_20261019.json
journal_flux_axial_mobility_tables_20261019.tex
journal_flux_axial_mobility_20261019.tex
journal_campaign_20261019.bib
journal_campaign_20261019_references.json
journal_flux_axial_mobility_20261019.pdf
```

The scientific suppliers import the hierarchy and local-flux model modules. The compact supplier additionally imports the canonical carrier and reads its result JSON; execution therefore follows the carrier, reusing the mass-three nucleon baseline while calculating the mass-six nucleon and compact six-quark branches at masses six and three. Its merger comparison is a finite canonical compact-channel result, not a hard/soft matching calculation. The carrier also reads the bytes of `rich_colour_spin_root_20261019.py` for its recorded source hash, so that file is included even though it is not imported or executed by reproduction. Its historical nine-link graph is not substituted for the paper's independent-loop bouquet. The hierarchy module's dormant standalone `run()` reads an older campaign input JSON; the supplied carrier/band/compact paths use its operators, never that standalone routine. Consequently that unrelated older input and its computation are intentionally not included or invoked. No additional external data files are required by the three supplied calculations. All selected scientific modules are flat beside the entrypoint, as required by their `HERE`-relative imports and source-hash reads.

The bibliography metadata file documents reference retrieval; it is not an executable input or an external service requirement. Historical source-provenance labels in the byte-preserved scientific sources/results are opaque identifiers, not credentials, capability tokens or required external services. All computation is local. Scientific files and recorded hashes are not silently redacted or transformed. Campaign narratives, private memory plans, session histories, instructions, mail, credentials, caches and unrelated research outputs are excluded. The builder rejects detected private absolute path roots in selected payload bytes rather than rewriting them. The numerical source/result provenance includes only the selected scientific files; no repository-wide history or operational dossier is bundled.
