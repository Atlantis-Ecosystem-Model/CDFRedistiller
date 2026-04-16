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
  
## Forcing File Timing on Restart

When restarting Atlantis from a checkpoint, the model continues with `bm->t`
set to the restart time (e.g. `bm->t = 63072000` seconds for a restart at
day 730). Whether each forcing source produces correct values on restart
depends on **how that forcing is indexed in the C code**.

### Restart-Safe Forcing (no action needed)

These forcing sources use `tsEval(ts, var_id, bm->t)` from the Atlantis time
series library, which interpolates by absolute time. Restart works correctly
as long as the `.ts` file covers the restart period:

| Forcing | Source file pattern | Code location |
|---|---|---|
| Forced effort | `Effortts*.data` in force.prm | `atManage.c:489-494` |
| Imposed catch | catch `.ts` files | `atHarvestImposedCatch.c:282` |
| Imposed discards | discard `.ts` files | `atHarvestDiscards.c:188` |
| Forced recruitment | `tsRecruit` `.ts` files | `atdemography.c:1247, 4587` |
| Environmental recruitment scalars | `Recruitment_enviro_forcing` | `atdemography.c:1686` |
| KWSR forcing | `KWSR_forcing` | `atdemography.c:1716` |
| Growth rate / FSPB / size scalers | various `.ts` files | `atecologyts.c`, `atannualbiology.c` |
| Linear mortality scaler | `tslinearMort` | `atq10.c:79` |
| pCO2 forcing | `tspCO2` `.ts` file | `atbiophysics.c:761` |
| Fuel cost forcing | fuel `.ts` files | `ateconindicator.c:285, 626` |
| Air temperature index | thermal index `.ts` | `atannualbiology.c:1266` |
| MPA closures (timeseries) | MPA `.ts` files | `atManageMPATS.c:265` |

For restart to work correctly with these, ensure your `.ts` files have time
entries covering the restart day onward. Files using `days since YYYY-MM-DD`
or `seconds since YYYY-MM-DD` as their time units are interpreted as absolute
calendar times. Atlantis automatically reconciles unit differences between
`.ts` files and the model via `tsNewTimeUnits`.

### Restart-Unsafe Forcing (currently breaks on restart)

The following forcing sources read records sequentially from netCDF files
using a step counter (`nextrec`) rather than absolute time. On restart,
`nextrec` initializes to 0 — so the model reads forcing values from the
**beginning of the file** while the simulation thinks it is at the restart
time. This causes a multi-year offset between the forcing applied and the
forcing intended.

| Forcing | Source | Code location |
|---|---|---|
| Hydrodynamic exchange | hydro `.nc` files | `athydromod.c:216, 396` |
| Temperature / salinity | tempsalt `.nc` files | `attempsalt.c:680, 820` |
| Ice | ice `.nc` files | `aticeIO.c:1067, 1184` |

The pattern in each case is:

```c
bm->hd.nextrec++;       // increment each timestep
...
bm->hd.nextrec = 0;     // initialize to 0 when file opens
```

There is no logic to seek forward to `bm->t` when the file is opened.

### Two Options for Fixing netCDF Forcing on Restart

**Option 1 — Modify `open_hydro` (and equivalents) to seek to `bm->t`** *(preferred, not yet implemented)*

After opening the netCDF file, read the `t` variable, find the index
where `t[i] >= bm->t`, and set `nextrec = i`. This makes the change once
and all future restarts work automatically without per-restart file
preparation. The same pattern applies to:
- `open_hydro` in `atphysics/athydromod.c`
- `open_phyprop` in `atphysics/attempsalt.c` (handles temperature, salinity, pH, wind, vertical mixing, light, noise, and tracer forcing)
- `Ice_Read_Time_Series` in `atphysics/aticeIO.c`

Estimated changes: ~10 lines of C per location, three locations.

**Option 2 — Pre-trim netCDF forcing files to start at the restart time** *(no C changes required)*

Use `ncks` to extract only the timesteps from the restart day onward:

```bash
# Trim hydro file (time in seconds since model epoch)
ncks -d t,63072000.,3153600000. original_hydro.nc trimmed_hydro.nc

# Trim tempsalt file
ncks -d t,63072000.,3153600000. original_tempsalt.nc trimmed_tempsalt.nc
```

Then point your restart `force.prm` at the trimmed files. The model still
starts `nextrec = 0`, but record 0 of the trimmed file now corresponds to
the restart time. This option requires maintaining a separate set of
forcing files per restart point.

### Recommendation

Option 1 is the long-term fix and removes any per-restart file
preparation. It has not yet been implemented. Until then, Option 2 (file
trimming) is the workaround for runs that include hydrodynamic, temperature,
salinity, or ice forcing on restart.


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