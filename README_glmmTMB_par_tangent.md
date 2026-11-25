# Parallel glmmTMB fitting (side note)

Notes on `R/glmmTMB_par.R` and the changes made to the glmmTMB candidate-fit scripts
(`model_fitting/HOR_cand_mods_run.R`, `TCJ_cand_mods_run.R`, `HOR_TCJ_cand_mods_run.R`).
Written 2026-11-25. Kept separate from the top-level `README.md` on purpose.

## What was added

`R/glmmTMB_par.R` holds four documented functions (roxygen comments). Each of the three
fitting functions starts its own worker cluster (12 workers by default, `n_workers` argument)
unless one is passed in `cl`.

| Function | Replaces | Used in |
|---|---|---|
| `dredge_par_glmmTMB()` | `MuMIn::dredge()` on a glmmTMB global model | HOR, TCJ |
| `get_models_par()` | `MuMIn::get.models(<dredge table>, subset = NA)` | HOR, TCJ |
| `fit_glmmTMB_par()` | the `lapply()` over formulas (`fit_HOR_TCJ_models()`) | HOR_TCJ |
| `glmm_cluster()` | starts the worker cluster used by the three above | (internal) |

In the scripts the original lines are commented out; directly below each is a
`source("R/glmmTMB_par.R")` (once per script, above the first new call) and the call to the new
function. Nothing else in the scripts was touched.

Implementation notes:

- Workers are started with `--no-init-file`. Started in the project folder they would each run
  `.Rprofile` and activate renv, which took about 20 s for 12 workers (now about 6 s including
  loading glmmTMB). The main session's library paths are passed to the workers explicitly.
- The data frame a model was fitted to is copied to the workers under the same name; the
  models are refitted there from their calls (`clusterApplyLB`, so a free worker takes the next
  model).
- `dredge_par_glmmTMB()` uses `MuMIn::dredge(cluster = )`. MuMIn numbers models from 0 in a
  cluster run and from 1 in a serial run, so the function shifts the row names (and the names
  of `model.calls`) to the serial numbering. The scripts pick models out by row name.
  `MuMIn::dredge()` also evaluates the data variable from inside its own namespace when it
  runs on a cluster, so the data frame must be a global variable (it is, in these scripts).

## Timing

Elapsed seconds for one `source()` of each script (12 workers, 16 logical cores):

| Script | Serial (HEAD) | Parallel refits only | Parallel dredge and refits |
|---|---|---|---|
| `HOR_cand_mods_run.R` | 366 | 237 | 103 |
| `TCJ_cand_mods_run.R` | 186 | 137 | 68 |
| `HOR_TCJ_cand_mods_run.R` | 31 | 18 | not applicable (no dredge) |

The parallel dredge accounted for roughly half of the saving in the HOR and TCJ scripts, which
is why `dredge_par_glmmTMB()` was kept. The gain for HOR_TCJ is small because the cluster
start-up takes a large share of its run time.

## Correctness check

The parallel results were compared with the stored serial results in `output/` (before the
`fixed` change below): AIC and fixed-effect estimates differ by 0, the model calls are identical,
and the AIC tables are equal with identical row names (HOR, TCJ and HOR_TCJ).

## The `fixed` term in the dredge calls (changed)

The calls used `fixed = c("barrier.exit.hor.fac")` (and `.tcj.fac`). glmmTMB names its terms with
a component prefix, so MuMIn did not find the term, warned that it does not exist in the global
model, and fixed nothing. The dredge therefore ran over 70 models, 16 of them without the
barrier. The four live calls now use `fixed = c("cond(barrier.exit.hor.fac)")` and
`fixed = c("cond(barrier.exit.tcj.fac)")`: 54 models per set, all containing the barrier main
effect. The commented-out original calls were left as they were.

Rerun into a test folder (HOR 98 s, TCJ 61 s) and compared with the stored results:

- **HOR:** the dAIC-filtered sets (6 models, unscaled and standardized) contain the same models
  with the same AIC, dAIC and k. Only the row numbers change, because the model numbering
  changes (for example 10, 28, 26, 74, 12, 14 become 5, 14, 13, 37, 6, 7).
- **TCJ:** the set changes. Unscaled and standardized sets agree with each other.

| | Before (`fixed` not applied) | After (barrier fixed) |
|---|---|---|
| Models in the dredge | 70 | 54 |
| Models in the dAIC-filtered set | 9 | 10 |
| Best model | CVP + log.VNS (AIC 2129.90, no barrier) | CVP + log.VNS + SWP + barrier + SWP:barrier (AIC 2130.54) |
| Models without a barrier term in the set | 5 | 0 |
| Shared with the other set | 4 | 4 |

The 4 shared models are the barrier models of the old set; their AIC values are unchanged. The
5 old models dropped from the set are all barrier-free (CVP + log.VNS, CVP, CVP + flength,
CVP + flength + log.VNS, CVP + log.VNS + SWP). 6 models enter that are new to the set. dAIC and the
AIC weights are relative to the best model, so they shift when the best model changes (dAIC of
the shared models rises by about 0.64).

Anything downstream of the TCJ dAIC-filtered set (`TCJ_d2_mods.rds` and what reads it) will change
when the pipeline is run with the new `fixed` term. The stored `output/` files still hold the
earlier results.
