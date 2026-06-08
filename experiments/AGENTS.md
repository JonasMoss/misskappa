# Experiments

Rules for one-off research / engineering investigations under
`experiments/`. The goal is that "write an experiment about ..." prompts
produce one focused, reproducible folder instead of a pile of notes.

## What an experiment is

A folder that answers **one** research or engineering question about
misskappa. It can have a large design grid, multiple estimators, or many
simulation cells, but those pieces must serve the same question. If the
question changes, start a new folder.

Experiments are advisory. They can motivate library work, paper
content, or new tests, but they are not themselves parity gates and do
not block releases. When an experiment is paper-worthy, it migrates into
the appropriate `papers/<slug>/` folder and the experiment folder stays
as a reproducibility record.

Anything still in development that isn't yet a question to be answered
belongs in `dev/notes/` instead.

## Lifecycle

The canonical experiment tree is organized by longevity promise:

- `studies/` — reusable or paper-facing studies that users or papers may
  inspect later.
- `workbench/` — active work in progress.
- `probes/` — diagnostics, stress tests, pilots, and validation checks.
- `archive/pre-redesign/` — frozen records whose runners are not maintained
  against the current R API.

Root-level `experiments/NN-slug/` directories are compatibility landing pages
only. New references should use the canonical lifecycle path listed in
`experiments/INDEX.md`.

## Directory shape

Start new work under `experiments/workbench/` with the next numeric prefix and
a short kebab-case slug. Promote it to `studies/`, `probes/`, or
`archive/pre-redesign/` once its status is clear, and update
`experiments/INDEX.md` in the same change:

```text
experiments/
  workbench/
    NN-topic-slug/
      report.qmd                # the narrative; the readable artifact
      run_experiment.R          # the compute entry point
      .gitignore                # ignores Quarto output + local resources
      results/
        .gitignore              # keeps results/ present; ignores raw runs
      R/                        # optional experiment-local helpers
      scripts/                  # optional secondary runners for the same question
      resources/                # optional local-only inputs; always ignored
```

Canonical experiment folders use the same internal shape:

```text
experiments/<lifecycle>/
  NN-topic-slug/
    report.qmd                # the narrative; the readable artifact
    run_experiment.R          # the compute entry point
    .gitignore                # ignores Quarto output + local resources
    results/
      .gitignore              # keeps results/ present; ignores raw runs
    R/                        # optional experiment-local helpers
    scripts/                  # optional secondary runners for the same question
    resources/                # optional local-only inputs; always ignored
```

Required files: `report.qmd`, `run_experiment.R`, `.gitignore`, and
`results/.gitignore`. Do not add a per-experiment `README.md` — the
report is the readable document and `Rscript run_experiment.R --help`
is the command reference. The exception is root-level legacy stubs, which
contain only a `README.md` pointing to the canonical path.

## Runners

`run_experiment.R` does the expensive work. Rendering the report must be
fast: read result files, reshape modest summaries, make small plots,
stop with a clear message when results are missing. Do not simulate,
fit, install packages, or regenerate fixtures from `report.qmd`.

Runners should:

- support `--help`;
- support a cheap smoke path (e.g., `--smoke`, `--reps 1`);
- accept a deterministic seed (`--seed-base` for simulations);
- create `results/` as needed;
- write `results/metadata.csv` with command arguments, seed base,
  session/package versions, and the design dimensions needed to
  interpret the report;
- write rectangular CSV outputs with stable column names;
- print the paths written at the end.

Runners load the misskappa R package directly (`library(misskappa)`).
If shared helpers accumulate, put them in
`experiments/_support/` as a small installable R package; do not put
estimator / inference logic there.

## Reports

Reports are professional and concise. Default shape:

1. Question
2. Short Answer
3. Evidence
4. Caveats
5. Reproduce

Lead with the point. Keep method detail only when it is needed to trust
the result. Reports read from `results/`; if files are absent, stop with
a command the reader can run.

Quarto defaults:

```yaml
format:
  html:
    toc: true
    number-sections: false
execute:
  echo: false
  warning: false
  message: false
```

## Tables and figures

Tables in reports are presentation objects, not data dumps.

- Human-readable headers (`Median Time (ms)`, not `median_elapsed_ms`).
- Short enough to fit on one screen. Summarise; do not paste the full
  design grid.
- Long tables live in `results/*.csv`, mentioned in the prose.
- Round numbers for interpretation; keep raw precision in the CSV.
- Prefer one clear plot over several near-duplicates.

For style of tables / figures / math / prose in *any* misskappa write-up,
follow `papers/ipw/STYLE.md` (cloned into each `papers/<slug>/STYLE.md`)
to the extent it fits a report.

## Style

Methods-developer voice: direct, specific, modest about claims. Say
what the experiment can and cannot establish. Separate statistical
conclusions from engineering diagnostics.

When a finding changes implementation state or the active backlog,
update `dev/notes/todo.md` or open a new entry there.
