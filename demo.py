# Self-contained demo: isotope envelopes by peptide size, and charge detection at two resolutions.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from isotopes import (isotope_envelope, monoisotopic_mass, mz, simulate_spectrum,
                      pick_peaks, charge_from_spacing)

AA = list("ACDEFGHIKLMNPQRSTVWY")


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)
    rng = np.random.default_rng(3)
    rows = []
    for length in [6, 10, 15, 20, 30, 50, 80, 120, 200, 400]:
        for rep in range(5):
            seq = "".join(rng.choice(AA, length))
            env = isotope_envelope(seq, peaks=60)
            rows.append({"length": length, "mono_mass": monoisotopic_mass(seq),
                         "mono_fraction": env[0], "most_abundant_offset": int(np.argmax(env))})
    sizes = pd.DataFrame(rows).groupby("length").mean().round(3).reset_index()

    pep = "".join(rng.choice(AA, 22))
    checks = []
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.8))
    for length, c in zip([10, 50, 200], ["#1f77b4", "#ff7f0e", "#c0504d"]):
        env = isotope_envelope("".join(rng.choice(AA, length)), peaks=14)
        axes[0].bar(np.arange(14) + {10: -0.27, 50: 0, 200: 0.27}[length], env, 0.27, color=c, label=f"{length} residues")
    axes[0].set_xlabel("isotope peak (M + n)")
    axes[0].set_ylabel("relative abundance")
    axes[0].set_title("The monoisotopic peak shrinks with size")
    axes[0].legend()

    for ax, res in zip(axes[1:], [60_000, 1_500]):
        for z, c in [(2, "#2a7f62"), (3, "#9467bd")]:
            centre = mz(monoisotopic_mass(pep), z)
            grid = np.linspace(centre - 1, centre + 4, 6000)
            y = simulate_spectrum(pep, z, res, grid, rng=rng)
            peaks = pick_peaks(grid, y)
            z_hat = charge_from_spacing(peaks)
            checks.append({"resolution": res, "true_charge": z, "peaks_found": len(peaks), "charge_called": z_hat})
            ax.plot(grid - centre, y + (1.2 if z == 3 else 0), color=c, lw=1, label=f"z={z}: called z={z_hat}")
        ax.set_xlabel("m/z offset from monoisotopic peak")
        ax.set_yticks([])
        ax.set_title(f"Resolution {res:,}")
        ax.legend(fontsize=8)
    table = pd.DataFrame(checks)
    sizes.to_csv("results/summary.csv", index=False)
    table.to_csv("results/charge_calls.csv", index=False)
    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    print(sizes.to_string(index=False))
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
