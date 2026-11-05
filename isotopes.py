# Peptide isotope patterns: elemental composition, isotope envelope, m/z, and charge from peak spacing.
import numpy as np

# Residue formulas (C, H, N, O, S) without water.
RESIDUES = {
    "A": (3, 5, 1, 1, 0), "R": (6, 12, 4, 1, 0), "N": (4, 6, 2, 2, 0), "D": (4, 5, 1, 3, 0),
    "C": (3, 5, 1, 1, 1), "E": (5, 7, 1, 3, 0), "Q": (5, 8, 2, 2, 0), "G": (2, 3, 1, 1, 0),
    "H": (6, 7, 3, 1, 0), "I": (6, 11, 1, 1, 0), "L": (6, 11, 1, 1, 0), "K": (6, 12, 2, 1, 0),
    "M": (5, 9, 1, 1, 1), "F": (9, 9, 1, 1, 0), "P": (5, 7, 1, 1, 0), "S": (3, 5, 1, 2, 0),
    "T": (4, 7, 1, 2, 0), "W": (11, 10, 2, 1, 0), "Y": (9, 9, 1, 2, 0), "V": (5, 9, 1, 1, 0),
}
MONO = {"C": 12.0, "H": 1.00782503, "N": 14.0030740, "O": 15.9949146, "S": 31.9720707}
# Natural abundances by nominal mass offset (+0, +1, +2, ...).
ISOTOPES = {"C": [0.9893, 0.0107], "H": [0.999885, 0.000115], "N": [0.99636, 0.00364],
            "O": [0.99757, 0.00038, 0.00205], "S": [0.9499, 0.0075, 0.0425, 0.0, 0.0001]}
PROTON = 1.00727646688
NEUTRON_SPACING = 1.00335  # average spacing between isotope peaks of peptides


def composition(seq):
    c, h, n, o, s = (sum(RESIDUES[a][k] for a in seq) for k in range(5))
    return {"C": c, "H": h + 2, "N": n, "O": o + 1, "S": s}


def monoisotopic_mass(seq):
    return sum(MONO[e] * k for e, k in composition(seq).items())


def isotope_envelope(seq, peaks=12):
    """Relative abundance of M, M+1, M+2, ... by repeated convolution of element patterns."""
    env = np.array([1.0])
    for element, count in composition(seq).items():
        pattern = np.array(ISOTOPES[element])
        for _ in range(count):
            env = np.convolve(env, pattern)[:peaks]
    return env / env.sum()


def mz(mass, charge):
    return (mass + charge * PROTON) / charge


def simulate_spectrum(seq, charge, resolution, mz_grid, noise=0.002, rng=None):
    rng = rng or np.random.default_rng()
    env = isotope_envelope(seq)
    centre0 = mz(monoisotopic_mass(seq), charge)
    fwhm = centre0 / resolution
    sigma = fwhm / 2.3548
    y = np.zeros_like(mz_grid)
    for k, a in enumerate(env):
        y += a * np.exp(-0.5 * ((mz_grid - centre0 - k * NEUTRON_SPACING / charge) / sigma) ** 2)
    return y / y.max() + rng.normal(0, noise, len(mz_grid))


def pick_peaks(mz_grid, y, min_height=0.05, min_separation=0.15):
    """Smoothed local maxima at least min_separation apart, centroided by intensity."""
    step = mz_grid[1] - mz_grid[0]
    w = max(int(min_separation / step), 2)
    box = max(w // 4, 1)
    smooth = np.convolve(y, np.ones(box) / box, mode="same")
    cents = []
    for i in range(w, len(y) - w):
        if smooth[i] > min_height and smooth[i] == smooth[i - w:i + w + 1].max():
            seg = slice(i - w // 2, i + w // 2 + 1)
            weights = np.clip(smooth[seg], 0, None)
            cents.append((np.sum(mz_grid[seg] * weights) / weights.sum(), smooth[i]))
    return np.array(cents).reshape(-1, 2)


def charge_from_spacing(centroids):
    if len(centroids) < 2:
        return None
    spacing = np.median(np.diff(np.sort(centroids[:, 0])))
    return int(round(NEUTRON_SPACING / spacing))
