"""Generate the paper's figures, tables and provenance from executed model results.

Run with the scientific Python environment; no supplier computation is repeated.
All continuous curves below are declared analytical retained-order expressions.
"""
from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

HERE = Path(__file__).resolve().parent
STEM = "journal_flux_axial_mobility"
DATE = "20261019"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table_macro(name: str, caption: str, label: str, columns: str,
                header: str, rows: list[str]) -> str:
    return ("\\newcommand{\\" + name + "}{%\n\\begin{table}[tbp]\n"
            "\\centering\\small\n\\caption{" + caption + "}\\label{" + label + "}\n"
            "\\begin{tabular}{" + columns + "}\\toprule\n" + header + "\\\\\n\\midrule\n"
            + "\\\\\n".join(rows) + "\\\\\n\\bottomrule\\end{tabular}\n\\end{table}\n}\n")


def main() -> None:
    carrier_path = HERE / "rich_journal_flux_carrier_20261019_result.json"
    band_path = HERE / "rich_journal_baryon_band_20261019_result.json"
    compact_path = HERE / "rich_journal_compact_hierarchy_20261019_result.json"
    carrier = json.loads(carrier_path.read_text())
    band = json.loads(band_path.read_text())
    compact = json.loads(compact_path.read_text())
    compact_rows = compact["comparisons"]
    canonical = sorted((r for r in carrier["rows"]
                        if r["loop_representation"] == "PW155_degree2" and r["I"] == .5),
                       key=lambda r: -r["mass_m"])
    delta_rows = {r["mass_m"]: r for r in carrier["rows"]
                  if r["loop_representation"] == "PW155_degree2" and r["I"] == 1.5}
    weak = next(r for r in canonical if r["mass_m"] == 3.)
    response = weak["inelastic_axial"]
    reference = response["matched_25_state_reference"]
    exterior = weak["canonical_degree3_exterior"]
    compact_weak = next(r for r in compact_rows if r["mass_m"] == 3.)
    compact_state = compact_weak["compact_N6"]
    nucleon_radius = float(np.sqrt(exterior["ground_state_H_exterior_norm2"]
                                   + weak["unprojected_H_eigen_residual"]**2))
    compact_radius = float(np.sqrt(compact_state["canonical_degree3_ground_exterior"]["H_exterior_norm2"]
                                   + compact_state["unprojected_H_eigen_residual"]**2))
    nucleon_radius_up = float(np.ceil(nucleon_radius*1e6)/1e6)
    compact_radius_up = float(np.ceil(compact_radius*1e6)/1e6)
    compact_gap_down = float(np.floor(compact_weak["compact_minus_separated"]*1e7)/1e7)
    rotor_inclusion = {
        "mass_m": 3.0,
        "meaning": "Numerical evaluation of sector-restricted full-local-rotor spectral inclusion; existence of compatible Gamma6 nucleon and I0 T1 six-quark eigenlevels, not global minima or adiabatic branch identification.",
        "formula": "dist(E,spectrum(h_sector)) <= sqrt(interior_residual^2 + exact_exterior_norm2)",
        "nucleon_center": weak["energy"], "compact_center": compact_state["energy"],
        "nucleon_interior_residual": weak["unprojected_H_eigen_residual"],
        "compact_interior_residual": compact_state["unprojected_H_eigen_residual"],
        "nucleon_exterior_norm2": exterior["ground_state_H_exterior_norm2"],
        "compact_exterior_norm2": compact_state["canonical_degree3_ground_exterior"]["H_exterior_norm2"],
        "nucleon_radius": nucleon_radius, "compact_radius": compact_radius,
        "nucleon_interval": [weak["energy"]-nucleon_radius, weak["energy"]+nucleon_radius],
        "compact_interval": [compact_state["energy"]-compact_radius, compact_state["energy"]+compact_radius],
        "hierarchy_margin": compact_weak["compact_minus_separated"]-compact_radius-2*nucleon_radius,
        "nucleon_radius_rounded_up": nucleon_radius_up,
        "compact_radius_rounded_up": compact_radius_up,
        "compact_gap_rounded_down": compact_gap_down,
        "rounded_hierarchy_margin": compact_gap_down-compact_radius_up-2*nucleon_radius_up,
        "numerical_status": "Floating-point residual norms with conservative decimal presentation; not directed interval arithmetic or a formally certified enclosure."}
    mass_rows = band["matched_mass_curve"]
    c_n = band["propagation"]["dimensionless_C_Omega2"]
    c_d = band["Delta_plus_propagation"]["dimensionless_C_Omega2"]
    eta = np.linspace(0, max(r["t_over_Delta"] for r in mass_rows), 250)
    x = eta**2 / 48
    ka = np.linspace(-np.pi, np.pi, 301)
    band_curves = [{"m": r["m"], "ka": ka.tolist(),
                    "gap_g2": (r["Delta_N_gap_k"]["constant"] +
                               r["Delta_N_gap_k"]["cos_ka_coefficient"]*np.cos(ka)).tolist()}
                   for r in mass_rows]
    moment_keys = ["inelastic_strength_mu0", "energy_weighted_mu1", "energy_squared_mu2"]
    excess = [100*(response[k]/reference[k]-1) for k in moment_keys]
    channels = exterior["exterior_response_by_field_channel"]
    channel_keys = ["three_loops_fundamental", "two_loops_degree2_plus1",
                    "single_loop_10_bar10", "single_loop_15_bar15"]
    provenance_names = {carrier_path.name, band_path.name, compact_path.name, Path(__file__).name,
                        "rich_journal_baryon_band_20261019.py"}
    provenance_names.update(carrier["source_sha256"])
    provenance_names.update(compact["source_sha256"])
    provenance = {name: digest(HERE/name) for name in sorted(provenance_names)}
    data = {"schema": 1, "description": "Source-bound scientific figure and table data",
            "state_semantics": "Source-connected Gamma6/Gamma8 three-quark and T1 six-quark carrier eigenbranches, not proved global ground states. Original supplier ground_* keys denote retained-branch or trial-state diagnostics.",
            "operator_semantics": "Canonical three- and six-quark rows and response/exterior moments use isolated local bouquets with spatial hopping suppressed: h=H_loc and h2=P2*h*P2. Compact/separated comparisons retain independent P2 regulators per cell. Transport rows use the spatial bouquet-tree hopping expansion.",
            "source_sha256": provenance,
            "carrier_recorded_source_sha256": carrier["source_sha256"],
            "compact_recorded_source_sha256": compact["source_sha256"],
            "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                        "matplotlib": matplotlib.__version__},
            "canonical_nucleon_rows": canonical, "canonical_delta_rows": list(delta_rows.values()),
            "transport_rows": mass_rows,
            "compact_hierarchy_rows": compact_rows,
            "full_rotor_spectral_inclusion": rotor_inclusion,
            "analytical_curves": {"meaning": "Retained-order predictions, not new diagonalizations",
                "definitions": {"x": "eta^2/48", "kappa_N_over_bare": "1+C_N*Omega^2*x",
                    "kappa_D_over_bare": "1+C_D*Omega^2*x",
                    "gap": "stored constant + stored cos(ka) coefficient*cos(ka)",
                    "two_layer_axial": "10/9+5*Z0/9"},
                "eta": eta.tolist(), "x": x.tolist(),
                "kappa_N_over_bare": (1+c_n*x).tolist(),
                "kappa_D_over_bare": (1+c_d*x).tolist(),
                "curvature_contrast": ((c_d-c_n)*x).tolist(), "band_gaps": band_curves},
            "matched_response_percent_excess": dict(zip(moment_keys, excess)),
            "exterior_channel_percent": {k:100*channels[k]/sum(channels.values()) for k in channel_keys}}
    (HERE/f"{STEM}_data_{DATE}.json").write_text(json.dumps(data, indent=2)+"\n")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False, "pdf.fonttype": 42,
                         "figure.constrained_layout.use": True})
    with PdfPages(HERE/f"{STEM}_figures_{DATE}.pdf",
                  metadata={"Title": "Chromomagnetic flux, mobility and axial response", "Author": "Anonymous"}) as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.1))
        ax = axes[0]
        ax.plot(eta, 1+c_n*x, label=r"$\kappa_N/\kappa_0$")
        ax.plot(eta, 1+c_d*x, label=r"$\kappa_\Delta/\kappa_0$")
        for r in mass_rows:
            bare = 12*r["t_NR"]**3*r["a"]**2 / band["root"]["Delta"]**2
            for label, color in [("N", "C0"), ("Delta", "C1")]:
                ax.scatter(r["t_over_Delta"], r["bands"][label]["curvature"]/bare, color=color, s=24)
        ax.set(xlabel=r"$\eta=t/\Delta$", ylabel="Curvature / bare curvature", title="(a) Opposite flux-induced mobility shifts")
        ax.legend(frameon=False)
        ax = axes[1]
        for curve in band_curves:
            ax.plot(np.array(curve["ka"])/np.pi, curve["gap_g2"], label=f'$m={curve["m"]:g}$')
        ax.axhline(0, color=".5", lw=.7)
        ax.set(xlabel=r"$ka/\pi$", ylabel=r"$E_\Delta(k)-E_N(k)$", title="(b) Retained-order momentum gap")
        ax.legend(frameon=False)
        pdf.savefig(fig); plt.close(fig)
        fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.1))
        z = np.linspace(.2, 1, 200)
        axes[0].plot(z, 10/9+5*z/9, color=".5", label="Vacuum / one-source law")
        for r in canonical:
            axes[0].scatter(r["vacuum_seed_overlap_Z0"], r["axial_neutral_tau3_sigma3"], s=38)
            axes[0].annotate(f'$m={r["mass_m"]:g}$', (r["vacuum_seed_overlap_Z0"], r["axial_neutral_tau3_sigma3"]), xytext=(4,-13), textcoords="offset points")
        axes[0].set(xlabel=r"Vacuum weight $Z_0$", ylabel=r"$g_A$", title="(a) Canonical two-source axial vertex")
        axes[0].legend(frameon=False, loc="upper left")
        axes[1].bar([r"$\mu_0$", r"$\mu_1$", r"$\mu_{2,P_2}$"], excess, color=["C0", "C1", "C2"])
        for j, value in enumerate(excess): axes[1].text(j, value+.35, f"{value:.2f}%", ha="center")
        axes[1].set(ylabel="Excess over E- and gA-matched reference (%)", ylim=(0,23), title=r"(b) Inelastic discrimination at $m=3$")
        pdf.savefig(fig); plt.close(fig)
        fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.1))
        labels = ["Three fundamental loops", "Degree-two + fundamental", r"Single loop $10,\overline{10}$", r"Single loop $15,\overline{15}$"]
        percentages = [data["exterior_channel_percent"][k] for k in channel_keys]
        axes[0].barh(labels[::-1], percentages[::-1], color="C3")
        axes[0].set(xlabel="Share of exterior centered-response norm (%)", title=r"(a) Missing degree-three channels, $m=3$")
        axes[1].bar(["Degree 0", "Degree 1", "Degree 2"], exterior["ground_probability_by_total_representation_degree"], color="C0")
        axes[1].set(yscale="log", ylabel="Carrier-state probability", title="(b) Retained representation-degree weights")
        pdf.savefig(fig); plt.close(fig)
    tables = "% Generated from component JSONs; do not hand-edit numerical entries.\n"
    rows = [f'{r["m"]:g} & {r["t_over_Delta"]:.6f} & {r["bands"]["N"]["curvature"]:.9f} & {r["bands"]["Delta"]["curvature"]:.9f} & {r["Delta_N_gap_k"]["minimum"]:.9f}' for r in mass_rows]
    tables += table_macro("TransportTable", "Transport predictions at retained magnetic and hopping orders. All rows have $a=1$, $k_E=0.27$. The final column uses the leading magnetic rest energy.", "tab:transport", "rrrrr", r"$m$ & $t/\Delta$ & $\kappa_N$ & $\kappa_\Delta$ & $\min_k(E_\Delta-E_N)$", rows)
    rows = [f'{r["mass_m"]:g} & {r["energy"]:.9f} & {delta_rows[r["mass_m"]]["energy"]-r["energy"]:.9f} & {r["axial_neutral_tau3_sigma3"]:.9f} & {r["vacuum_seed_overlap_Z0"]:.9f} & {r["probability_law_defect"]:.9f}' for r in canonical]
    tables += table_macro("CanonicalTable", "Source-connected carrier levels and vertices of the canonical total-degree-two Galerkin Hamiltonian, with $a=1$, $k_E=0.27$. Energies are in the model's common energy units.", "tab:canonical", "rrrrrr", r"$m$ & $E_N$ & $E_\Delta-E_N$ & $g_A$ & $Z_0$ & $D_A$", rows)
    labels_mom = [r"$\mu_0$", r"$\mu_1$", r"$\mu_{2,P_2}$", r"$\mu_1/\mu_0$"]
    keys_mom = moment_keys + ["centroid_mu1_over_mu0"]
    rows = [f'{label} & {response[k]:.9f} & {reference[k]:.9f} & {100*(response[k]/reference[k]-1):+.3f}\\%' for label,k in zip(labels_mom,keys_mom)]
    tables += table_macro("ResponseTable", "Local neutral inelastic axial response at $m=3$ and $t=0$. The reference is fitted to both the canonical carrier energy and elastic axial vertex; both final-isospin sectors are included. The second moment is internal to $P_2$.", "tab:response", "lrrr", r"Observable & Canonical $P_2$ & Matched reference & Difference", rows)
    texlabels = ["Three fundamental loops", "Degree-two plus fundamental loop", r"Single $10,\overline{10}$ loop", r"Single $15,\overline{15}$ loop"]
    rows = [f'{label} & {channels[k]:.10f} & {data["exterior_channel_percent"][k]:.3f}\\%' for label,k in zip(texlabels,channel_keys)]
    tables += table_macro("ExteriorTable", "Exact exterior action of the untruncated local $h$ on the retained $m=3$ trial state at $t=0$: positive contributions to $\\|Q_2h(A-g_A)\\Psi_2\\|^2$. These are not spectral moments of an uncomputed full ground state.", "tab:exterior", "lrr", "Channel & Norm squared & Share", rows)
    rows = [f'{r["mass_m"]:g} & {r["compact_N6"]["energy"]:.10f} & {r["separated_pn_I0_T1_energy"]:.10f} & {r["compact_minus_separated"]:.10f} & {r["leading_compact_minus_separated"]:.10f}' for r in compact_rows]
    tables += table_macro("CompactTable", "Local compact six-quark channel versus two independent nucleons in the matched $I=0$, $T_1$, charge-one channel. Both occupied cells use the same canonical $P_2$ regulator; $a=1$, $k_E=0.27$. The leading column is $(28/3)g^2/\\Omega$, not a fitted penalty.", "tab:compact", "rrrrr", r"$m$ & $E_6$ & $2E_N$ & $E_6-2E_N$ & Leading", rows)
    for name, value, digits in [
            ("RotorNRadius", nucleon_radius, 12), ("RotorSixRadius", compact_radius, 12),
            ("RotorNRadiusUp", nucleon_radius_up, 6), ("RotorSixRadiusUp", compact_radius_up, 6),
            ("RotorCompactGapDown", compact_gap_down, 7),
            ("RotorNumericalMargin", rotor_inclusion["hierarchy_margin"], 10),
            ("RotorConservativeMargin", rotor_inclusion["rounded_hierarchy_margin"], 7)]:
        tables += "\\newcommand{\\" + name + "}{" + f"{value:.{digits}f}" + "}\n"
    (HERE/f"{STEM}_tables_{DATE}.tex").write_text(tables)
    print(json.dumps({"figure_pages": 3, "canonical_rows": len(canonical),
                      "transport_rows": len(mass_rows), "compact_rows": len(compact_rows),
                      "spectral_inclusion_margin": rotor_inclusion["rounded_hierarchy_margin"],
                      "source_sha256": provenance}, indent=2))


if __name__ == "__main__":
    main()
