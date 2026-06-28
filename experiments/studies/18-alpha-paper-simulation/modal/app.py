"""Modal harness for the alpha three-sim (study 18).

Fans ONE container per design cell (categories, dgp, mechanism, n). Each runs
run_experiment.R for that cell (B reps) and writes its own summary/replicates to
a shared Volume. The per-fit time + memory guards SIGKILL any stalled or
ballooning cat-FIML fit and record it (killed / killed_mem) rather than crashing
the cell; Modal retries any cell that dies wholesale; seeds are stored per
replicate so a retry is deterministic. combine.R concatenates all cells into the
final summary.csv + replicates.csv.

Usage (from this directory, with the Modal CLI authenticated):

    modal run app.py                              # full run, B=1000, 90 cells
    modal run app.py --reps 20 --run-id smoke     # quick smoke
    modal run app.py --reps 200 --run-id pilot

Pull results:

    modal volume get alpha18-results <run-id>/summary.csv ./<run-id>-summary.csv
    modal volume get alpha18-results <run-id>/replicates.csv ./<run-id>-replicates.csv
"""

from pathlib import Path

import modal

if modal.is_local():
    HERE = Path(__file__).resolve()
    EXP = HERE.parents[1]                 # 18-alpha-paper-simulation/
    REPO_ROOT = HERE.parents[4]           # misskappa/
    RPKG = REPO_ROOT / "r-package"
else:
    EXP = Path("/exp")
    RPKG = Path("/build/r-package")

# rocker/r-ver ships a pinned R with the Posit binary repo, so Rcpp/RcppEigen
# install as binaries and the C++17 toolchain is present. misskappa is compiled
# from source in-image so cat-FIML / nt_fiml_cs are live (and carry the
# Louis-memory fix). Exclude host build artifacts so the in-image gcc recompiles.
image = (
    modal.Image.from_registry("rocker/r-ver:4.4.1", add_python="3.11")
    .run_commands("Rscript -e 'install.packages(c(\"Rcpp\", \"RcppEigen\"))'")
    .add_local_dir(str(RPKG), "/build/r-package", copy=True,
                   ignore=["**/*.o", "**/*.so"])
    .run_commands("R CMD INSTALL /build/r-package")
    .add_local_file(str(EXP / "run_experiment.R"), "/exp/run_experiment.R")
    .add_local_file(str(EXP / "modal" / "combine.R"), "/exp/combine.R")
)

app = modal.App("alpha18")
vol = modal.Volume.from_name("alpha18-results", create_if_missing=True)

METHODS = "pairwise,cat_em,normal_fiml,nt_fiml_cs,listwise"
MECHS = ["mcar30", "mar_anchor30"]
NGRID = [100, 250, 500, 1000, 2000]
TIMEOUT_S = "30"        # per-fit wall-clock kill
MEM_MB = "3000"         # per-fit RSS kill (Louis mem is bounded; this is insurance)

# Two category configs: cats=5 covers Sims A (4-item) + B (6-item) DGPs;
# cats=7 covers Sim C (the 7-point challenge, 6-item DGPs).
CONFIGS = [
    (5, ["ordinal_tau4", "ordinal_congeneric4", "ordinal_twofactor4",
         "ordinal_tau6", "ordinal_congeneric6", "ordinal_twofactor6"]),
    (7, ["ordinal_tau6", "ordinal_congeneric6", "ordinal_twofactor6"]),
]


def cells():
    """Enumerate (idx, cats, dgp, mech, n). idx drives a unique per-cell seed."""
    idx = 0
    for cats, dgps in CONFIGS:
        for dgp in dgps:
            for mech in MECHS:
                for n in NGRID:
                    idx += 1
                    yield (idx, cats, dgp, mech, n)


@app.function(image=image, volumes={"/vol": vol}, timeout=3 * 60 * 60,
              cpu=2.0, memory=4096, retries=3)
def run_cell(idx: int, cats: int, dgp: str, mech: str, n: int,
             run_id: str, reps: int) -> int:
    """One design cell: B reps at a single (categories, dgp, mechanism, n)."""
    import subprocess

    out_dir = f"/vol/{run_id}/cells/c{cats}_{dgp}_{mech}_n{n}"
    cmd = [
        "Rscript", "/exp/run_experiment.R",
        "--reps", str(reps), "--n-grid", str(n), "--categories", str(cats),
        "--dgps", dgp, "--mechanisms", mech, "--methods", METHODS,
        "--fit-timeout", TIMEOUT_S, "--fit-memory-mb", MEM_MB,
        "--seed-base", str(idx * 10_000_000), "--out-dir", out_dir,
    ]
    print("RUN", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    vol.commit()
    return idx


@app.function(image=image, volumes={"/vol": vol}, timeout=60 * 60,
              cpu=2.0, memory=8192)
def combine(run_id: str) -> str:
    """Concatenate every cell's summary/replicates into the final CSVs."""
    import subprocess

    vol.reload()
    out_dir = f"/vol/{run_id}"
    subprocess.run(["Rscript", "/exp/combine.R", "--out-dir", out_dir], check=True)
    vol.commit()
    return out_dir


@app.local_entrypoint()
def main(reps: int = 1000, run_id: str = "full"):
    args = [(i, c, d, m, n, run_id, reps) for (i, c, d, m, n) in cells()]
    print(f"Launching {len(args)} cells, reps={reps}, run_id={run_id}")
    # Blocks until all cells finish (each committed its own results).
    done = list(run_cell.starmap(args))
    print(f"Cells complete: {len(done)}/{len(args)}")
    out_dir = combine.remote(run_id)
    print(f"Combined into volume alpha18-results at {out_dir}.")
    print("Pull results with:")
    print(f"  modal volume get alpha18-results {run_id}/summary.csv ./{run_id}-summary.csv")
    print(f"  modal volume get alpha18-results {run_id}/replicates.csv ./{run_id}-replicates.csv")
