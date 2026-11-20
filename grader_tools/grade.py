"""Standalone grader for the wind-blown bubble task (implements grader.md).

Usage:
    python grader_tools/grade.py <ref_dir> [<ref_dir> ...]

Writes <ref_dir>/grades/ containing the per-task grade files and grade_summary.md.
Outputs are graded strictly against the format requested in prompt.md; a file that
does not match it counts as missing (grade 0) and the reason is listed in the summary.
"""
import glob
import os
import re
import sys

import numpy as np
from astropy import units as u
from astropy import constants as ac

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ad_wind_spec import AdiabaticWind  # noqa: E402

# Task parameters (from prompt.md)
SAM_EXPERIMENTS = {"F100": 100.0, "F10": 10.0, "F1": 1.0}  # name -> n_H [cm^-3]
MDOT = 1e-4 * u.Msun / u.yr
LWIND = 1e38 * u.erg / u.s
SAM_RFB = 100 * u.AU
SAM_T_KYR = 100.0
SIM_NH = 100.0
SIM_RFB = 0.5 * u.pc
N_RESOLUTIONS = 3
N_SNAPSHOTS = 5
CONV_VARIABLES = ["temperature", "density", "velocity", "pressure"]
CONV_EXPECTED_INDEX = {"fo": -1.0, "so": -2.0}
ORDERS = ["fo", "so"]
CSV_HEADER = ["r", "rho", "v", "P", "T"]
GRADE_MAX = 6.0


class FormatError(Exception):
    """An output file exists but does not follow the requested format."""


def log_grade(x):
    """min(max(-log10 x, 0), 6). A missing value (None or nan) scores 0."""
    if x is None or not np.isfinite(x):
        return 0.0
    if x == 0:
        return GRADE_MAX
    return float(min(max(-np.log10(x), 0.0), GRADE_MAX))


def make_model(n_H, rfb):
    rho0 = (1.4 * ac.m_p * n_H / u.cm**3).to(u.g / u.cm**3)
    return AdiabaticWind(Lwind=LWIND, Mdotw=MDOT, rfb=rfb.to(u.pc), rho0=rho0)


# ---------------------------------------------------------------------------
# Strict readers
# ---------------------------------------------------------------------------
def read_profile_csv(path):
    """Comma-separated; the first line is exactly CSV_HEADER (prompt.md, Output File Formats)."""
    with open(path) as f:
        header = f.readline().strip()
    if header != ",".join(CSV_HEADER):
        raise FormatError(f"{path}: first line is '{header[:80]}', expected '{','.join(CSV_HEADER)}'")
    try:
        data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2, comments=None)
    except ValueError as e:
        raise FormatError(f"{path}: could not parse data rows ({e})")
    if data.shape[1] != len(CSV_HEADER):
        raise FormatError(f"{path}: {data.shape[1]} data columns, expected {len(CSV_HEADER)}")
    return {name: data[:, i] for i, name in enumerate(CSV_HEADER)}


def read_number_lines(path, n_per_line):
    """Lines of exactly n_per_line space-separated numbers; no headers or comments.
    Blank lines (e.g. a trailing newline) are ignored."""
    rows = []
    for k, line in enumerate(open(path), start=1):
        if not line.strip():
            continue
        tokens = line.split()
        try:
            values = [float(t) for t in tokens]
        except ValueError:
            raise FormatError(f"{path}: line {k} is not numeric: '{line.strip()[:80]}'")
        if len(values) != n_per_line:
            raise FormatError(f"{path}: line {k} has {len(values)} numbers, expected {n_per_line}")
        rows.append(values)
    if not rows:
        raise FormatError(f"{path}: file is empty")
    return np.array(rows)


def read_two_numbers(path):
    """A single line: outer radius [pc], then P/k_B [K cm^-3]."""
    rows = read_number_lines(path, 2)
    if len(rows) != 1:
        raise FormatError(f"{path}: expected a single line, found {len(rows)}")
    return rows[0]


def read_times(path):
    """One time [kyr] per line; line k is snapshot k."""
    return list(read_number_lines(path, 1)[:, 0])


def read_convergence(path):
    """One line per resolution: N, L1 (first order), L1 (second order)."""
    return read_number_lines(path, 3)


# ---------------------------------------------------------------------------
# Comparisons
# ---------------------------------------------------------------------------
def profile_l1(model, prof, t_kyr):
    """Mean relative error of rho, v, P over 0 < r < R_bubble(t), skipping q_ref = 0."""
    t = t_kyr * u.kyr
    r = prof["r"]
    inside = (r > 0) & (r < model.radius(t).to_value(u.pc))
    r_in = r[inside] * u.pc
    reference = {
        "rho": model.density_profile(r_in, t).to_value(u.g / u.cm**3),
        "v": model.velocity_profile(r_in, t).to_value(u.km / u.s),
        "P": model.pressure_profile(r_in, t).to_value(u.K / u.cm**3),
    }
    l1 = {}
    for key, q_ref in reference.items():
        q = prof[key][inside]
        ok = q_ref != 0
        l1[key] = float(np.mean(np.abs(q[ok] - q_ref[ok]) / np.abs(q_ref[ok]))) if ok.any() else np.nan
    return l1


def grade_sam(ref_dir, notes):
    """Returns {experiment: {test: grade}}; writes <name>_sam_grade.txt."""
    tests = ["radius", "pressure", "rho", "v", "P"]
    labels = ["radius(t) vs .txt", "pressure(t) vs .txt", "density_profile",
              "velocity_profile", "pressure_profile"]
    out_dir = os.path.join(ref_dir, "sam_outputs")
    results = {}
    for name, n_H in SAM_EXPERIMENTS.items():
        l1 = dict.fromkeys(tests, np.nan)
        if os.path.isdir(out_dir):
            model = make_model(n_H, SAM_RFB)
            t = SAM_T_KYR * u.kyr
            try:
                R, P = read_two_numbers(os.path.join(out_dir, f"{name}.txt"))
                R_ref = model.radius(t).to_value(u.pc)
                P_ref = model.pressure(t).to_value(u.K / u.cm**3)
                l1["radius"], l1["pressure"] = abs(R - R_ref) / R_ref, abs(P - P_ref) / P_ref
            except (FormatError, OSError) as e:
                notes.append(f"SAM {name}.txt: {e}")
            try:
                l1.update(profile_l1(model, read_profile_csv(os.path.join(out_dir, f"{name}.csv")), SAM_T_KYR))
            except (FormatError, OSError) as e:
                notes.append(f"SAM {name}.csv: {e}")
        grades = {k: log_grade(v) for k, v in l1.items()}
        results[name] = grades
        with open(os.path.join(ref_dir, "grades", f"{name}_sam_grade.txt"), "w") as f:
            f.write(f"# {name}: -log10(relative-error L1 norm), clamped to [0, {GRADE_MAX:g}]\n")
            f.write("# " + " | ".join(labels) + "\n")
            f.write(" ".join(f"{grades[k]:.4f}" for k in tests) + "\n")
            f.write("# raw L1: " + " ".join(f"{l1[k]:.4e}" for k in tests) + "\n")
    if not os.path.isdir(out_dir):
        notes.append("sam_outputs/ missing: all SAM grades are 0")
    return results


def grade_sim(ref_dir, notes):
    """Returns {order: [row, ...]}; writes sim_grade_<order>.txt. Missing snapshots and
    resolutions (up to N_RESOLUTIONS x N_SNAPSHOTS per order) are listed as 0."""
    sim_dir = os.path.join(ref_dir, "sim_outputs")
    model = make_model(SIM_NH, SIM_RFB)
    pattern = re.compile(r"F100_sim_(fo|so)_N(\d{5})_snap(\d{3})\.csv$")
    found = {}  # (order, N) -> {snap: path}
    for path in glob.glob(os.path.join(sim_dir, "F100_sim_*_snap*.csv")):
        m = pattern.search(os.path.basename(path))
        if m and int(m[3]) >= 1:
            found.setdefault((m[1], int(m[2])), {})[int(m[3])] = path
        else:
            notes.append(f"sim: ignored {os.path.basename(path)} (name must be "
                         "F100_sim_VV_NXXXXX_snapYYY.csv, zero-padded, snapshots from 001)")
    resolutions = sorted({N for _, N in found})
    if not os.path.isdir(sim_dir):
        notes.append("sim_outputs/ missing: all simulation grades are 0")

    times = {}
    for N in resolutions:
        try:
            times[N] = read_times(os.path.join(sim_dir, f"F100_sim_N{N:05d}_times.txt"))
        except (FormatError, OSError) as e:
            notes.append(f"sim N={N} times: {e}")

    results = {}
    for order in ORDERS:
        rows = []
        for N in resolutions:
            snaps = found.get((order, N), {})
            graded = 0
            for snap in sorted(snaps):
                t_kyr, l1 = np.nan, dict.fromkeys(["rho", "P", "v"], np.nan)
                idx = snap - 1  # line k of the times file is snapshot k
                if N in times and idx < len(times[N]):
                    t_kyr = times[N][idx]
                    if t_kyr <= 0:
                        continue  # initial condition: not graded, not missing
                    try:
                        l1 = profile_l1(model, read_profile_csv(snaps[snap]), t_kyr)
                    except (FormatError, OSError) as e:
                        notes.append(f"sim {order} N={N} snap {snap:03d}: {e}")
                elif N in times:
                    notes.append(f"sim {order} N={N} snap {snap:03d}: no matching line in times file")
                rows.append((N, f"{snap:03d}", t_kyr, l1))
                graded += 1
            if graded < N_SNAPSHOTS:
                notes.append(f"sim {order} N={N}: {N_SNAPSHOTS - graded} of {N_SNAPSHOTS} snapshots missing")
            for _ in range(max(N_SNAPSHOTS - graded, 0)):
                rows.append((N, "missing", np.nan, dict.fromkeys(["rho", "P", "v"], np.nan)))
        if len(resolutions) < N_RESOLUTIONS:
            notes.append(f"sim {order}: {N_RESOLUTIONS - len(resolutions)} of {N_RESOLUTIONS} resolutions missing")
        for _ in range(max(N_RESOLUTIONS - len(resolutions), 0)):
            for _ in range(N_SNAPSHOTS):
                rows.append((0, "missing", np.nan, dict.fromkeys(["rho", "P", "v"], np.nan)))
        results[order] = rows
        with open(os.path.join(ref_dir, "grades", f"sim_grade_{order}.txt"), "w") as f:
            f.write(f"# F100 simulation ({order}) vs AdiabaticWind(rfb={SIM_RFB}): "
                    f"-log10(mean relative-error L1 over 0 < r < R_bubble(t)), clamped to [0, {GRADE_MAX:g}]\n")
            f.write("# N  snap  t_kyr  density  pressure  velocity | raw L1: density  pressure  velocity\n")
            for N, snap, t_kyr, l1 in rows:
                g = [log_grade(l1[k]) for k in ["rho", "P", "v"]]
                f.write(f"{N:6d} {snap:>7s} {t_kyr:8.3f}  {g[0]:.4f} {g[1]:.4f} {g[2]:.4f} | "
                        f"{l1['rho']:.4e} {l1['P']:.4e} {l1['v']:.4e}\n")
    return results


def grade_convergence(ref_dir, notes):
    """Returns a 2x4 array of grades (rows fo, so; columns CONV_VARIABLES)."""
    conv_dir = os.path.join(ref_dir, "convergence")
    indices = np.full((len(ORDERS), len(CONV_VARIABLES)), np.nan)
    if not os.path.isdir(conv_dir):
        notes.append("convergence/ missing: all convergence grades are 0")
    else:
        for j, var in enumerate(CONV_VARIABLES):
            try:
                data = read_convergence(os.path.join(conv_dir, f"{var}_L1_norm.txt"))
            except (FormatError, OSError) as e:
                notes.append(f"convergence {var}: {e}")
                continue
            for i in range(len(ORDERS)):
                N, l1 = data[:, 0], data[:, i + 1]
                if len(N) < 2 or np.any(N <= 0) or np.any(l1 <= 0):
                    notes.append(f"convergence {var} {ORDERS[i]}: need >= 2 positive points to fit")
                    continue
                indices[i, j] = np.polyfit(np.log10(N), np.log10(l1), 1)[0]
    expected = np.array([[CONV_EXPECTED_INDEX[o]] for o in ORDERS])
    grades = np.vectorize(log_grade)(np.abs(indices - expected))
    with open(os.path.join(ref_dir, "grades", "convergence_grades.txt"), "w") as f:
        f.write(f"# -log10|fitted power-law index - expected|, clamped to [0, {GRADE_MAX:g}]; "
                "expected -1 (first order), -2 (second order)\n")
        f.write("# rows: first order, second order; columns: " + " ".join(CONV_VARIABLES) + "\n")
        for row in grades:
            f.write(" ".join(f"{g:.4f}" for g in row) + "\n")
        f.write("# fitted indices (nan = missing):\n")
        for row in indices:
            f.write("# " + " ".join(f"{x:+.4f}" for x in row) + "\n")
    return grades, indices


def read_code_quality(ref_dir):
    """Optional manual score: first number in grades/code_quality.txt, else None."""
    path = os.path.join(ref_dir, "grades", "code_quality.txt")
    if not os.path.exists(path):
        return None
    for line in open(path):
        if line.strip() and not line.lstrip().startswith("#"):
            try:
                return float(line.replace(",", " ").split()[0])
            except ValueError:
                return None
    return None


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
def write_summary(ref_dir, sam, sim, conv_grades, conv_indices, code_quality, notes):
    sam_total = sum(sum(g.values()) for g in sam.values())
    sim_totals = {o: sum(log_grade(l1[k]) for *_, l1 in rows for k in ["rho", "P", "v"])
                  for o, rows in sim.items()}
    conv_total = float(conv_grades.sum())
    final = sam_total + sum(sim_totals.values()) + conv_total + (code_quality or 0.0)

    lines = [f"# Grade summary: `{os.path.basename(os.path.abspath(ref_dir))}`", "",
             f"Every grade is min(max(-log10 x, 0), {GRADE_MAX:g}) of a relative-error L1 norm "
             "(or of a power-law-index difference), so higher is better and missing outputs score 0.", "",
             "## Semi-analytic model", "",
             "| Experiment | radius | pressure | density profile | velocity profile | pressure profile | total |",
             "| --- | --- | --- | --- | --- | --- | --- |"]
    for name, g in sam.items():
        lines.append(f"| {name} | " + " | ".join(f"{v:.3f}" for v in g.values())
                     + f" | {sum(g.values()):.3f} |")
    lines += ["", f"SAM total: **{sam_total:.3f}** (max {GRADE_MAX * 5 * len(sam):g})", "",
              "## Simulation (vs `AdiabaticWind`, rfb = 0.5 pc)", "",
              "| Order | snapshots graded | missing | total |", "| --- | --- | --- | --- |"]
    for o, rows in sim.items():
        n_missing = sum(1 for r in rows if r[1] == "missing")
        lines.append(f"| {o} | {len(rows) - n_missing} | {n_missing} | {sim_totals[o]:.3f} |")
    lines += ["", "Per-snapshot grades are in `sim_grade_fo.txt` and `sim_grade_so.txt`.", "",
              "## Convergence", "",
              "| Order | " + " | ".join(CONV_VARIABLES) + " |",
              "| --- |" + " --- |" * len(CONV_VARIABLES)]
    for i, o in enumerate(ORDERS):
        cells = [f"{g:.3f} (index {x:+.2f})" if np.isfinite(x) else f"{g:.3f} (missing)"
                 for g, x in zip(conv_grades[i], conv_indices[i])]
        lines.append(f"| {o} | " + " | ".join(cells) + " |")
    lines += ["", f"Convergence total: **{conv_total:.3f}**", "", "## Code quality", "",
              f"Manual score from `grades/code_quality.txt`: **{code_quality:g}**" if code_quality is not None
              else "Not graded (no `grades/code_quality.txt`). See README.md for how to add it.",
              "", "## Final grade", "",
              f"SAM {sam_total:.3f} + sim fo {sim_totals['fo']:.3f} + sim so {sim_totals['so']:.3f} "
              f"+ convergence {conv_total:.3f}"
              + (f" + code quality {code_quality:g}" if code_quality is not None else ""), "",
              f"**Final grade: {final:.3f}**", ""]
    if notes:
        lines += ["## Format problems and missing outputs", ""] + [f"- {n}" for n in notes] + [""]
    with open(os.path.join(ref_dir, "grades", "grade_summary.md"), "w") as f:
        f.write("\n".join(lines))
    return final


def grade_directory(ref_dir):
    os.makedirs(os.path.join(ref_dir, "grades"), exist_ok=True)
    notes = []
    sam = grade_sam(ref_dir, notes)
    sim = grade_sim(ref_dir, notes)
    conv_grades, conv_indices = grade_convergence(ref_dir, notes)
    code_quality = read_code_quality(ref_dir)
    final = write_summary(ref_dir, sam, sim, conv_grades, conv_indices, code_quality, notes)
    print(f"{ref_dir}: final grade {final:.3f} ({len(notes)} format/missing notes) "
          f"-> {os.path.join(ref_dir, 'grades', 'grade_summary.md')}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for d in sys.argv[1:]:
        if not os.path.isdir(d):
            sys.exit(f"not a directory: {d}")
        grade_directory(d)
