# Style

How misskappa papers look and read. Companion to `paper/AGENTS.md` (which
defines manuscript structure and paper-specific direction). When you draft
text or build a table or figure, this file is the reference.

## Audience Profiles

Each paper's AGENTS.md declares a `Profile:` line. Three options:

| Profile | Target journals | Math | Theorems | Default body (pages) | Brief slot |
|---|---|---|---|---|---|
| **Math** | Psychometrika, BJMSP, JEBS | dense | yes, in body | **17–22** | Short Note, <10 pp |
| **Applied** | SEM, MBR, Psychological Methods (also EPM for measurement work) | moderate, in displays | rare | **10–14** | Notes & Comments, <10 pp |
| **Tool** | BRM | light, tutorial OK | almost never | **11–14** | Short Report, <10 pp |

Pages are total published pages (manuscript including refs), from a
2020–2025 OpenAlex pull of articles ≥ 5 pages across the listed journals.
Targets are the **lower quartile** of each journal cluster, on the
default-aim-short rule. Real medians: Psychometrika 26, BJMSP 25, JEBS 29
(the longest cluster), MBR 19, Psychological Methods 19, EPM 25, SEM 14,
BRM 18. A paper may overshoot when there is content to justify it; it
should not overshoot by default.

What each profile expects beyond length:

- **Math.** Theorem/Proof environments fine. Symbol density expected.
  Figures sparing; one or two carefully chosen plots. Discussion short.
  Prose closer to a clean mathematical exposition than to a tutorial.
- **Applied.** Theorems only as math displays plus paragraph statement,
  never as `Theorem 1` environments. Simulation study expected.
  Discussion in plain prose. Reproducibility against irrCAC is central.
- **Tool.** Tutorial framing OK. Software description with usage examples.
  Validation studies stand in for theory. Code and data effectively
  required at submission.

A `Brief slot` paper drops to roughly half the default length and
typically carries one of: a closed-form result, a parity correction, a
small software release, or a single-finding empirical note. The
manuscript should be submitted explicitly to the journal's brief
category, not as a short regular article.

## Tables

APA shape, `booktabs` always, no exceptions:

- **Rules.** `\toprule`, `\midrule`, `\bottomrule` only. Never `\hline`,
  never `|` between columns.
- **Caption above the tabular**, ending in a period, self-contained. Lead
  with the verb. Name sample size, design dimensions, and the headline
  finding when there is one. A reader who skims the captions should still
  follow the paper.
- **Spanner headers** via `\cmidrule(lr){a-b}` instead of repeating
  column-name prefixes. If three columns are all `Bias (SD)` for three
  estimators, put `Bias (SD)` once over a `\cmidrule`-spanned group and
  let the three sub-columns name the estimators.
- **Decimal alignment** via `siunitx` `S` columns. One number format per
  column, never per row. Round for interpretation; keep raw precision in
  the corresponding CSV under `results/`.
- **Headers in human Title Case** (`Median Time (ms)`), never raw
  identifiers (`median_elapsed_ms`).
- **At most ~12 columns.** Above that, transpose or split. Spanner groups
  count once.
- **Group repeated row labels** with an extra `\midrule` and a group label;
  do not repeat the same string down the leftmost column.
- **Notes below the tabular** in a `\parbox` matching table width.
  General Note first, one short sentence on what is shown. Symbol
  footnotes second, alphabetic superscripts in order of appearance. Notes
  describe *what*, not *how*; method details belong in the body text.

Worked example:

```latex
\begin{table}[t]
\caption{Bias and Monte Carlo SD of Conger's $\kappa_C$ under DGP B
(non-exchangeable raters, MCAR missingness). IPW corrects the bias that
AC and Gwet inherit from rater-specific observation rates.}
\label{tab:sim-conger}
\centering
\begin{tabular}{l S[table-format=2.4] S[table-format=2.4] S[table-format=2.4]}
\toprule
        & \multicolumn{3}{c}{Bias (SD), $n = 600$} \\
\cmidrule(lr){2-4}
Estimator & {Conger} & {Fleiss} & {Brennan--Prediger} \\
\midrule
AC   & 0.1039 & 0.1041 & 0.1058 \\
IPW  & 0.0027 & 0.0030 & 0.0046 \\
Gwet & 0.1046 & 0.1048 & 0.1058 \\
\bottomrule
\end{tabular}
\par\smallskip
\footnotesize\textit{Note.} 500 replications per cell.
\end{table}
```

What to refuse:

- Vertical rules and `\hline`.
- Footnote bushes longer than the table itself.
- A "results" column with three rows in different units.
- Captions that say only `Table 3: Timings.`
- A 30-row design grid pasted whole; summarise with top-N or grouped
  slices and refer the reader to the CSV.

## Figures

APA shape, monochrome-safe, paper-consistent:

- **One theme** across every figure in every paper. `theme_classic()` base
  with the modifications listed below.
- **Fixed palette per dimension**, reused across figures and across papers.
  Estimators are always the same colour wherever they appear (AC = black,
  IPW = blue, Gwet = red, FIML = green; agreement weighting families get
  their own stable assignment). Use the Okabe-Ito palette below as the
  default; document the per-dimension assignment once at the top of the
  figure script.
- **Linetypes vary alongside colour.** A black-and-white print of any
  figure must still be readable. The "headline" series is solid; comparators
  are dashed, dotted, dot-dash.
- **Captions** follow the table rules: above-the-figure in APA, but in
  practice LaTeX places them below; lead with the verb, self-contained,
  headline finding inline when there is one.
- **Vector formats only.** PDF for the manuscript, no PNG for plots.
- **Standard widths.** Single-column `3.5"` (~89 mm), two-column `7"`
  (~178 mm). Default to single-column. Heights chosen for a roughly 3:2
  or 4:3 aspect ratio, never wider than 7" or taller than 5".
- **No chartjunk.** No 3D, no shadows, no gradients, no chart background
  shading, no per-bar texture. Gridlines only where they aid reading
  (use sparing horizontal gridlines for value comparison; never both
  axes).
- **Axes** labelled with units, sensible breaks, no superfluous tick
  marks. Numbers respect the table number-format rules.
- **Legends** inside the plot if there is room and they do not overlap
  data; otherwise to the right. Title-Case labels, never raw identifiers.

The Okabe-Ito palette (colourblind-safe, eight colours):

```
black   #000000
orange  #E69F00
sky     #56B4E9
green   #009E73
yellow  #F0E442
blue    #0072B2
red     #D55E00
purple  #CC79A7
```

A typical per-dimension assignment for one paper:

```r
paper_palette <- list(
  estimator = c(ac = "#000000", ipw = "#0072B2",
                gwet = "#D55E00", fiml = "#009E73"),
  weight    = c(identity = "#000000", linear = "#0072B2",
                quadratic = "#D55E00", ordinal = "#009E73")
)
paper_linetype <- list(
  estimator = c(ac = "solid", ipw = "dashed",
                gwet = "dotted", fiml = "dotdash")
)
```

Pick the assignment once per paper, write it at the top of the figure
script, and reuse.

## Citations

APA throughout. `natbib` with `plainnat` (or `apalike`) is the default
and matches what the existing manuscript uses. Stay with it.

- `\citet{key}` for textual citations: `Author (Year) showed...`
- `\citep{key}` for parenthetical: `...(Author, Year).`
- `\citep[][p.~123]{key}` for page references.
- `\citep[e.g.,][]{key}` for "e.g." prefixes.
- `\citep{a,b,c}` for multi-citation; let `apalike` sort them.

`.bib` hygiene:

- Author names in proper case (`Bollen, K. A.`), never all-caps. Use
  protective braces (`{B}ollen`) only where the bibstyle would otherwise
  lowercase wrongly.
- Journal names full or consistently abbreviated within a paper. Pick
  one style for the whole `.bib` and stick to it.
- Include DOIs when available; URLs only for materials without a DOI.
- Distinguish preprint and published versions explicitly; do not silently
  cite the preprint when the journal version exists.
- One `.bib` per paper folder. Do not share `.bib` files across papers;
  copy the entries you need.

## Math

A short notation contract that holds across all misskappa papers unless
the per-paper AGENTS.md overrides for a documented reason:

- Italic lowercase for scalars and vectors: `$x$`, `$\theta$`, `$y$`.
- Italic uppercase for matrices: `$X$`, `$\Sigma$`, `$\Lambda$`.
- Bold reserved for explicit collections when the distinction matters:
  `$\boldsymbol{\theta}$` for the full parameter vector when individual
  `$\theta_j$` also appear and need disambiguation. Do not bold by default.
- Multi-letter operators via `\operatorname{}`: `\operatorname{vec}`,
  `\operatorname{tr}`, `\operatorname{diag}`, `\operatorname*{arg\,min}`.
  Never write `vec(X)` as upright text; always as an operator.
- Equation numbers only on equations the prose references. Unnumbered
  displays for transitional algebra.
- Do not display a single-symbol expression. `$x$` stays inline.
- Use `\,` for thin spaces in integrals and product names (`\,dx`,
  `\arg\,min`).

**Symbol glossary.** Every paper carries a short Symbol Glossary section
near the front of the manuscript, listing each symbol once with its
meaning. Prevents the "p means both p-value and number of indicators"
drift. Update the glossary whenever notation changes.

Theorems by profile:

- **Math.** `Theorem`, `Lemma`, `Proposition`, `Corollary`, `Proof`
  environments fine. Proofs in the body when short, in a short appendix
  when long. State assumptions explicitly above each statement.
- **Applied.** No `Theorem` environments. State the result as a numbered
  display plus a short paragraph naming assumptions and consequence.
  Citations carry the rigour.
- **Tool.** No formal mathematical statements unless they are immediate
  consequences of definitions.

## Prose

The voice we use, the words we do not.

### Voice

Pick "we" or impersonal per paper and stick with it. Recommended:

- **Math, Tool.** Impersonal: "The estimator solves...", "The
  simulation reports..."
- **Applied.** "We" is fine: "We estimate kappa via...", "We compare to
  irrCAC."

Past tense for what was done (methods, results). Present tense for what
holds (definitions, theorems, ongoing implications). Do not slip into
future tense ("we will show that") outside a paragraph that genuinely
forecasts later material.

### Structural rules

These hold across every profile:

- **No em-dashes.** Use a comma, a pair of parentheses, or two sentences.
  En-dashes in compound names (Golub-Pereyra) and numeric ranges are
  correct and stay.
- **No semicolons.** Split into two sentences, or use a comma.
- **Few colons.** A colon may introduce a genuine list. Do not use one
  for a dramatic pause before a punchline.
- **No antithesis.** Sentences that negate a foil and then assert the
  real point ("It is not X, it is Y", "not one decision but several",
  "this is not a failure but a demarcation") almost never carry
  information the positive statement does not. State the point directly,
  or drop the sentence.
- **No short blogger assertions.** Punchy standalone sentences that
  announce or characterise the text instead of advancing it
  ("The recommendation is short.", "This paper is an engineering
  study.") belong in the trash. Open paragraphs with substance.
- **Methods-developer voice.** Direct, specific, modest about claims.
  Say what the work can and cannot establish. Separate statistical
  conclusions from engineering diagnostics.

### Words and phrases to refuse

LLM-isms and stock academic filler that should never appear in a
misskappa paper:

- **Empty meta-comments.** "It is worth noting that...", "It is
  important to note that...", "Notably,", "Importantly,",
  "Interestingly,", "It is interesting to note that...". The reader can
  decide what is interesting. Say the thing.
- **Filler verbs.** "delve into", "dive into", "explore", "leverage",
  "utilize", "facilitate", "demonstrate". Use "show", "use", "let".
- **Filler adjectives.** "robust", "powerful", "elegant", "intuitive",
  "comprehensive", "novel" when they describe the work being presented.
  Show the property; do not assert it.
- **Vague quantifiers.** "various", "numerous", "myriad", "plethora",
  "a wide range of". Name the number or the cases.
- **Empty connectives at sentence start.** "Moreover,", "Furthermore,",
  "Additionally,", "In addition,". Either the sentence connects on its
  own or it should be merged with the previous one.
- **Bookend phrases.** "In this section we will...", "Having discussed
  X, we now turn to Y", "In conclusion,". The section structure is its
  own signal.
- **Hedge stacks.** "may potentially", "could possibly", "appears to
  perhaps", "tends to generally". Pick one hedge or none.
- **Triadic padding.** "X, Y, and Z" structures repeated across many
  consecutive sentences. Vary the rhythm. If you have two items, write
  two.
- **Clause chains.** "X is critical for Y, ensuring Z and providing W."
  Break into two sentences or drop the trailing clauses.
- **False precision.** "approximately 78.34%", "roughly 1.245$\times$".
  Round when hedging; do not hedge a precise number.
- **"Notably" stand-ins.** "Of particular interest", "It bears
  mentioning", "Of note", "Remarkably". Same problem, slightly different
  costume.

Mostly these compound: a sentence with two of them is wrong twice. When
in doubt, delete the sentence and re-read; if nothing breaks, the
sentence was filler.

### What good looks like

- A clear claim per sentence, with the verb early.
- Numbers come from the macros file, never typed by hand.
- Caveats stated in the same paragraph as the claim, not deferred.
- Method paragraphs end with what the choice buys, not with what it is.

## Build Hygiene for Style

Style consistency that the pipeline must enforce, not the author:

- **Every cited number is generated.** Prose numbers come from a single
  `tables/<slug>_stats.tex` file written by the paper's table script.
  Caption numbers come from the same file. The manuscript carries `??`
  fallbacks for unbuilt macros so it compiles before the tables exist.
- **Figure and table file slugs match in-text `\label{}` slugs.**
  `figures/sim_conger_bias.pdf` paired with
  `\label{fig:sim-conger-bias}` (kebab-case in labels, snake-case in
  filenames is fine, but the stem must match).
- **Tables and figures prefixed with the paper slug.**
  `tables/kappa-missing_sim_summary.tex`, not `tables/summary.tex`.
  Prevents cross-paper collisions if another paper ever lands in the
  repo, and makes grep useful.
- **Symbol glossary updated when notation changes.** A grep for `$p$`
  should agree with what the glossary says `p` means.
- **One R script writes one rectangular CSV.** A script that produces
  three artifacts produces three CSVs. Tables and figures read CSVs,
  never call simulation code.

## Notes on This File

These conventions hold by default. A paper may override a rule when the
paper-specific AGENTS.md explains why. Conventions that get overridden
in two or more papers should migrate back here.

The published-page numbers in the audience-profile table come from a
2020-2025 OpenAlex query (`type:article`, >=5 pages, filtered by journal
ISSN).
