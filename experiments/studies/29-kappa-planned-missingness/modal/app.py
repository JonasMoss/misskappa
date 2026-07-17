"""Modal harness for Experiment 29 (kappa planned-missingness master sim).

Builds an image with R + the local `misskappa` source (C++ backend compiled
in-image, so the co-observation guard is live), fans the replicate shards across
containers writing per-cell checkpoints to a shared Volume, then runs combine.R
to reduce all shards into the final replicates.csv + summary.csv.

Usage (from this directory, with the Modal CLI authenticated):

    modal run app.py                         # targeted, 32 shards, reps from mode
    modal run app.py --mode targeted --shard-count 64
    modal run app.py --mode smoke --shard-count 4 --reps 5 --run-id smoke-modal
    modal run app.py --mode targeted --shard-count 2 --reps 20 --run-id dryrun

Then pull results locally:

    modal volume get kappa29-results <run-id>/summary.csv ./<run-id>-summary.csv
    modal volume get kappa29-results <run-id>/replicates.csv ./<run-id>-replicates.csv

and render the report against the downloaded directory:

    KAPPA29_RESULTS_DIR=results/<run-id> quarto render ../report.qmd
"""

from pathlib import Path

import modal

# Local paths feed the image definition; inside the container app.py is
# mounted at a shallow path where parents[4] does not exist, so guard with
# is_local() (the image is already built there and these are never used).
if modal.is_local():
    HERE = Path(__file__).resolve()
    EXP = HERE.parents[1]                  # 29-kappa-planned-missingness/
    REPO_ROOT = HERE.parents[4]            # misskappa/
    RPKG = REPO_ROOT / "r-package"
else:
    EXP = Path("/exp")
    RPKG = Path("/build/r-package")

# rocker/r-ver ships a pinned R with the Posit binary repo (RSPM) preconfigured,
# so Rcpp/RcppEigen install as binaries and the C++17 toolchain is already there.
image = (
    modal.Image.from_registry("rocker/r-ver:4.4.1", add_python="3.11")
    .run_commands("Rscript -e 'install.packages(c(\"Rcpp\", \"RcppEigen\"))'")
    # Build-time copy of the package source -> compile + install misskappa.
    # Exclude local build artifacts: stale .o/.so compiled by the host's gcc
    # would be linked as "up to date" inside the image and fail to load there
    # (e.g. undefined __cxa_call_terminate from a newer host libstdc++).
    .add_local_dir(str(RPKG), "/build/r-package", copy=True,
                   ignore=["**/*.o", "**/*.so"])
    .run_commands("R CMD INSTALL /build/r-package")
    # Runtime mounts of the experiment scripts (editing these does NOT rebuild
    # the image). Must come after all run_commands.
    .add_local_file(str(EXP / "run_experiment.R"), "/exp/run_experiment.R")
    .add_local_file(str(EXP / "summarize.R"), "/exp/summarize.R")
    .add_local_file(str(EXP / "combine.R"), "/exp/combine.R")
)

app = modal.App("kappa29")
vol = modal.Volume.from_name("kappa29-results", create_if_missing=True)


@app.function(image=image, volumes={"/vol": vol}, timeout=6 * 60 * 60,
              cpu=2.0, memory=8192)
def run_shard(shard_index: int, shard_count: int, run_id: str, mode: str,
              reps: int) -> int:
    """Run one replicate shard: every design cell, this shard's replicate subset."""
    import subprocess

    out_dir = f"/vol/{run_id}"
    cmd = [
        "Rscript", "/exp/run_experiment.R", f"--{mode}",
        "--shard-index", str(shard_index),
        "--shard-count", str(shard_count),
        "--out-dir", out_dir,
        "--resume",            # reuse this shard's checkpoints on retry
        "--checkpoints-only",  # final reduce is combine.R's job
        "--progress",
    ]
    if reps > 0:
        cmd += ["--reps", str(reps)]
    print("RUN", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    vol.commit()
    return shard_index


@app.function(image=image, volumes={"/vol": vol}, timeout=60 * 60,
              cpu=2.0, memory=8192)
def combine(run_id: str) -> str:
    """Merge every shard's per-cell checkpoints into replicates.csv + summary.csv."""
    import subprocess

    vol.reload()  # see checkpoints committed by the shard containers
    out_dir = f"/vol/{run_id}"
    subprocess.run(["Rscript", "/exp/combine.R", "--out-dir", out_dir], check=True)
    vol.commit()
    return out_dir


@app.local_entrypoint()
def main(mode: str = "targeted", shard_count: int = 32, reps: int = 0,
         run_id: str = ""):
    if not run_id:
        run_id = mode
    args = [(i, shard_count, run_id, mode, reps) for i in range(1, shard_count + 1)]
    print(f"Launching {shard_count} shard(s) for mode={mode} run_id={run_id} "
          f"reps={'mode-default' if reps == 0 else reps}")
    # Blocks until all shards finish (and have committed their checkpoints).
    done = list(run_shard.starmap(args))
    print(f"Shards complete: {sorted(done)}")
    out_dir = combine.remote(run_id)
    print(f"Combined into volume kappa29-results at {out_dir}.")
    print("Pull results with:")
    print(f"  modal volume get kappa29-results {run_id}/summary.csv ./{run_id}-summary.csv")
    print(f"  modal volume get kappa29-results {run_id}/replicates.csv ./{run_id}-replicates.csv")
