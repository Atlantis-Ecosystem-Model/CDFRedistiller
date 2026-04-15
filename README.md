# cdf_redist — Atlantis CDF Redistributor

Python replacement for cdfDistiller (Beth Fulton, CSIRO, 2006).

Creates a restart init.nc from an Atlantis output.nc by copying the original
init file and surgically overwriting the dynamic variables with values from
the output at a chosen timestep.

## Requirements

- Python 3.6+
- netCDF4 (`pip install netCDF4`)
- numpy

## Quick Start

```bash
python3 cdf_redist.py -i Out.nc -o Init_restart.nc --init-file InitialCond.nc
```

## Usage

```
python3 cdf_redist.py -i <output.nc> -o <new_init.nc> --init-file <original_init.nc> [options]

Required:
  -i, --input       Atlantis output NetCDF file
  -o, --output      New init NetCDF file to create
  --init-file       Original init file (required for attributes and structure)

Options:
  -t, --timestep    Timestep index to extract (default: -1 = last)
  -v, --verbose     Print per-variable details
  --no-flip         Skip layer reordering (for testing)
  --force           Overwrite output file if it exists
```

## Examples

```bash
# Extract last timestep (default)
python3 cdf_redist.py -i Out.nc -o Init_restart.nc --init-file InitialCond.nc

# Extract a specific timestep
python3 cdf_redist.py -i Out.nc -o Init_restart.nc --init-file InitialCond.nc -t 5

# Verbose output
python3 cdf_redist.py -i Out.nc -o Init_restart.nc --init-file InitialCond.nc -v
```

## How It Works

1. Copies the original init file byte-for-byte (preserving all variables,
   attributes, fill values, masks, and ordering)
2. Opens the copy in read-write mode
3. Overwrites only the dynamic variables with output values at the chosen timestep:
   - All (t,b,z) tracer variables (bmtype='tracer') — with layer flip
   - dz — with layer flip
   - All (t,b) epibenthos variables
   - Physical (t,b) variables: numlayers, topk, sed*
   - Time variable
4. Leaves everything else untouched (volume, water, nominal_dz, porosity,
   habitat vars, disabled groups, all attributes)

## Time Handling

The script writes the time value from the selected output timestep into the
restart init.nc. This is the time **in seconds** at which the first run ended.
For example, if the first run covered 2 years (730 days), the time value in the
restart init.nc will be `63072000.0` (730 × 86400).

This preserved time value is critical because Atlantis uses it to index
time-dependent forcing, harvest schedules, and spawning calendars from the
original model start (t=0). Resetting time to zero would desynchronize all
time-dependent inputs.

The `--keep-time` flag has been removed. Time is always preserved from the
output file — this is the only correct behavior for restart runs.

## Run Configuration for Restart

After creating the restart init.nc, update the run configuration:

```bash
# In RunAtlantis.sh — use the restart init file:
atlantisMerged -i Init_restart.nc 0 -o Out.nc ...

# In Run.prm — tstop must cover the FULL timeline from t=0:
tstop      3650 day   # Total model time from original t=0, not remaining time

# toutstart should be set to the restart time or later:
toutstart  730 day    # First output at or after restart point
```

## Layer Ordering

Atlantis uses opposite layer ordering in input vs output files:

```
Input:  [L0, L1, L2, L3, 0, 0, SED]   — data first, zeros, sediment last
Output: [0, 0, L0, L1, L2, L3, SED]   — zeros first, data, sediment last
```

The tool automatically flips (t,b,z) variables from output to input convention.

## Age Structure (popratio)

For best results when restarting, restore the within-age-class structure
(per Beth Fulton's recommendation):

1. In the original run, set `flagagecheck 1` in Run.prm
2. Extract `popratioStock` vectors from the age check diagnostics
3. In the restart run's biol.prm, set `readin_popratio 1` and provide the vectors

Without this, Atlantis assumes uniform age distribution within each cohort,
causing a short burn-in period. Even with popratio, some burn-in is expected.

## Required Atlantis C Code Changes

Starting from a non-zero time exposes bugs in the spawning date initialization.
Four changes are required in the Atlantis C source code. All changes are
marked with `/* JMK CDFRedistiller */` comments.

### 1. Fix run-year calculation (`atbiolsetup.c`)

In `Ecology_Setup_Reproduction`, the number of run years must be calculated
from the remaining time, not total time.

```c
/* Original: */
int numRunYears = (int)(ceil(bm->tstop / 86400) / 365.0 + 1);

/* Fixed (JMK CDFRedistiller): */
int numRunYears = (int)(ceil((bm->tstop - bm->t) / 86400) / 365.0 + 1);
```

### 2. Offset spawning dates (`atannualbiology.c`, ~line 460)

In the spawning date setup loop (inside the `bcohort` loop), add the restart
day offset so that spawning dates are absolute rather than relative to t=0.

```c
/* Original: */
EMBRYO[sp].Spawn_Now[cohort][qid] += (yr_scalar * 365.0);
EMBRYO[sp].StartDay[cohort][lid]  += (yr_scalar * 365.0);
EMBRYO[sp].EndDay[cohort][lid]    += (yr_scalar * 365.0);

/* Fixed (JMK CDFRedistiller): */
int start_day_offset = (int)(bm->t / 86400.0); 
EMBRYO[sp].Spawn_Now[cohort][qid] += (yr_scalar * 365.0) + start_day_offset;
EMBRYO[sp].StartDay[cohort][lid]  += (yr_scalar * 365.0) + start_day_offset;
EMBRYO[sp].EndDay[cohort][lid]    += (yr_scalar * 365.0) + start_day_offset;
```

### 3. Skip spawn date backfill on restart (`atannualbiology.c`, ~line 505)

The backfill logic creates "previous year" spawning entries by subtracting
multiples of 365 from the current spawning date. On a restart, the dates are
already absolute, and the backfill creates entries before the restart time
which then overwrite the valid dates. Disable backfill when restarting.

```c
/* After the while (temp_val > 364) loop, before if (countback_spawn): */

/* JMK CDFRedistiller - no backfill needed when restarting, dates are already absolute */
if (bm->t > 0.0)
    countback_spawn = 0;

// Create the inserted countback cases
if (countback_spawn) {
    ...
```

### 4. Skip recruit date backfill on restart (`atannualbiology.c`, ~line 540)

Same issue as spawn backfill, but for recruitment dates.

```c
/* After the while (temp_val_end > 364) loop, before if (countback_recruit): */

/* JMK CDFRedistiller - no backfill needed when restarting, dates are already absolute */
if (bm->t > 0.0)
    countback_recruit = 0;

if (countback_recruit) {
    ...
```

### Known Limitations

- Does not yet work for migratory species. Migration dates likely need similar
  offset treatment.

## Companion R Code

This R code compares the terminal-year catch and biomass between a full
(uninterrupted) run and a restarted run to verify that the checkpoint
restart produces equivalent results.

```R
library(tidyverse)

# --- Configuration ---
base <- "path_to_out_dir"  # Parent output directory

full_run    <- file.path(base, "run_dir")  # Full continuous run
restart_run <- file.path(base, "run_dir")  # Stopped + restarted run

# Labels for output (adjust to match your setup)
label_full    <- "full"
label_restart <- "restart"

# --- Read data ---
catch_full    <- read_delim(file.path(full_run, "OutCatch.txt"), delim = " ", show_col_types = FALSE)
catch_restart <- read_delim(file.path(restart_run, "OutCatch.txt"), delim = " ", show_col_types = FALSE)

biom_full    <- read_delim(file.path(full_run, "OutBiomIndx.txt"), delim = " ", show_col_types = FALSE)
biom_restart <- read_delim(file.path(restart_run, "OutBiomIndx.txt"), delim = " ", show_col_types = FALSE)

# --- Terminal year values ---
terminal_catch_full <- catch_full |>
  slice_tail(n = 1) |>
  pivot_longer(-Time, names_to = "species", values_to = "catch") |>
  mutate(run = label_full)

terminal_catch_restart <- catch_restart |>
  slice_tail(n = 1) |>
  pivot_longer(-Time, names_to = "species", values_to = "catch") |>
  mutate(run = label_restart)

terminal_biom_full <- biom_full |>
  slice_tail(n = 1) |>
  pivot_longer(-Time, names_to = "species", values_to = "biomass") |>
  mutate(run = label_full)

terminal_biom_restart <- biom_restart |>
  slice_tail(n = 1) |>
  pivot_longer(-Time, names_to = "species", values_to = "biomass") |>
  mutate(run = label_restart)

# --- Compare catch ---
catch_compare <- bind_rows(terminal_catch_full, terminal_catch_restart) |>
  select(-Time) |>
  pivot_wider(names_from = run, values_from = catch) |>
  mutate(
    diff = .data[[label_restart]] - .data[[label_full]],
    pct_diff = 100 * diff / .data[[label_full]]
  ) |>
  filter(.data[[label_full]] > 0 | .data[[label_restart]] > 0) |>
  arrange(desc(abs(pct_diff)))

# --- Compare biomass ---
biom_compare <- bind_rows(terminal_biom_full, terminal_biom_restart) |>
  select(-Time) |>
  pivot_wider(names_from = run, values_from = biomass) |>
  mutate(
    diff = .data[[label_restart]] - .data[[label_full]],
    pct_diff = 100 * diff / .data[[label_full]]
  ) |>
  filter(.data[[label_full]] > 0 | .data[[label_restart]] > 0) |>
  arrange(desc(abs(pct_diff)))

cat("=== CATCH COMPARISON (terminal year) ===\n")
print(catch_compare, n = 50)

cat("\n=== BIOMASS COMPARISON (terminal year) ===\n")
print(biom_compare, n = 50)
```

## Author

Jacob Kasper (Marine and Freshwater Research Institute, Iceland)

Based on cdfDistiller by Beth Fulton (CSIRO), 2006