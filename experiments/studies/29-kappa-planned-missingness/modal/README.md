# Modal harness — Experiment 29

Runs the planned-missingness master simulation on [Modal](https://modal.com),
fanning replicate shards across containers and reducing them with `combine.R`.

## Prerequisites

- `pip install modal` and `modal token new` (one-time auth).
- The local `r-package/` source is baked into the image at build time, so the
  installed `misskappa` includes the co-observation guard. Editing the package
  triggers an image rebuild; editing `run_experiment.R` / `summarize.R` /
  `combine.R` does **not** (they are runtime mounts).

## Run

```bash
cd experiments/studies/29-kappa-planned-missingness/modal

# Paper-grade targeted run: 5 DGPs x 9 mechanisms x n{10,20,40,100} x 2000 reps.
# Always --detach: without it a local network blip stops the app mid-run
# (recoverable via --resume + same --run-id, but avoidable).
modal run --detach app.py --mode targeted --shard-count 32

# Cheap dry run first (2 shards, 20 reps) to validate the volume + combine path.
modal run app.py --mode targeted --shard-count 2 --reps 20 --run-id dryrun
```

`--shard-count` is the number of replicate shards (containers). Size it to the
cell/rep cost — calibrate from the smoke/screen wall-clock. Each shard runs every
design cell over its replicate subset and writes
`/<run-id>/checkpoints/cell-*-shardMofN.csv` to the `kappa29-results` Volume;
`combine.R` then reduces all shards into `summary.csv` + `replicates.csv`.

## Collect + report

```bash
modal volume get kappa29-results targeted/summary.csv ./results/targeted/summary.csv
modal volume get kappa29-results targeted/replicates.csv ./results/targeted/replicates.csv
# (also pull truth.csv, dgps.csv, mechanisms.csv, metadata.csv for the report)

KAPPA29_RESULTS_DIR=results/targeted quarto render ../report.qmd
```

## Notes

- Sharding is by **replicate**: distinct shards hold disjoint rep ids, so merged
  checkpoints are duplicate-free. `combine.R` also drops byte-identical rows
  defensively.
- Fixed-rater cells on the disconnected designs (`block2`, `ring2`, `random2`)
  are expected to fail the identifiability guard and show up as failures; the
  exchangeable `Fleiss-counts` track is what survives there.
