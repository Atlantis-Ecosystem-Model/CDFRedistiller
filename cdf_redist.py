#!/usr/bin/env python3
"""
cdf_redist.py — Atlantis CDF Redistributor

Python replacement for cdfDistiller (Beth Fulton, CSIRO, 2006).

Creates a restart init.nc from an Atlantis output.nc by copying the original
init file and surgically overwriting the dynamic variables with values from
the output at a chosen timestep.

Conservative approach:
  - Starts from a byte-perfect copy of the original init file
  - Only overwrites variables that change during a run
  - Preserves all attributes, fill values, masks, and variable ordering
  - Leaves volume, water, porosity, and other geometry vars untouched

TIME HANDLING:
  By default the time variable is reset to 0 because current Atlantis builds
  (as of 2025_10) have a bug in Init_Spawning that causes crashes when
  starting from a non-zero time. Use --keep-time to preserve the original
  output time value if this is fixed in a future build.

  When using t=0 (default), your forcing and harvest files must be adjusted
  to start from the restart point. Alternatively, if your forcing uses
  absolute dates, you may need to offset them.

LAYER ORDERING:
  Input:  [L0, L1, L2, L3, 0, 0, SED]   — data first, zeros, sediment last
  Output: [0, 0, L0, L1, L2, L3, SED]   — zeros first, data, sediment last

Usage:
  python3 cdf_redist.py -i Out.nc -o Init_restart.nc --init-file InitialCond.nc
  python3 cdf_redist.py -i Out.nc -o Init_restart.nc --init-file InitialCond.nc -t 5
  python3 cdf_redist.py -i Out.nc -o Init_restart.nc --init-file InitialCond.nc --keep-time

Repository: https://github.com/Atlantis-Ecosystem-Model/cdfdistiller
Author:     Jacob Kasper (Marine and Freshwater Research Institute, Iceland)
Based on:   cdfDistiller by Beth Fulton (CSIRO), 2006
Date:       2026
"""

import argparse
import sys
import os
import shutil
import numpy as np

try:
    import netCDF4 as nc
except ImportError:
    print("ERROR: netCDF4 required. pip install netCDF4")
    sys.exit(1)


def flip_layers_output_to_input(data_2d, numlayers, wcnz):
    """Flip (nbox,nz) from output [0,0,L0..Ln,SED] to input [L0..Ln,0,0,SED]."""
    nbox, nz = data_2d.shape
    flipped = np.zeros_like(data_2d)
    for b in range(nbox):
        nl = int(numlayers[b])
        if nl == 0:
            flipped[b, :] = data_2d[b, :]
            continue
        wc = data_2d[b, :wcnz]
        n_empty = wcnz - nl
        flipped[b, :nl] = wc[n_empty:]
        flipped[b, nl:wcnz] = wc[:n_empty]
        flipped[b, wcnz:] = data_2d[b, wcnz:]
    return flipped


def get_numlayers(ds_out, timestep, wcnz):
    if "numlayers" in ds_out.variables:
        nl = ds_out.variables["numlayers"][timestep, :]
        if hasattr(nl, "filled"):
            nl = nl.filled(0)
        return nl.astype(int)
    dz = ds_out.variables["dz"][timestep, :, :]
    if hasattr(dz, "filled"):
        dz = dz.filled(0.0)
    return np.sum(dz[:, :wcnz] > 0, axis=1).astype(int)


def time_to_human(t_val, t_units):
    try:
        return str(nc.num2date(t_val, t_units))
    except Exception:
        pass
    if "seconds" in t_units.lower():
        return f"~{float(t_val) / (365.25 * 24 * 3600):.1f} years from model start"
    return None


def parse_args():
    p = argparse.ArgumentParser(
        description="Create an Atlantis restart init from output + original init.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i Out.nc -o Init_restart.nc --init-file Init.nc
  %(prog)s -i Out.nc -o Init_restart.nc --init-file Init.nc -t 5
  %(prog)s -i Out.nc -o Init_restart.nc --init-file Init.nc -v
  %(prog)s -i Out.nc -o Init_restart.nc --init-file Init.nc --keep-time
        """,
    )
    p.add_argument("-i", "--input", required=True, help="Atlantis output.nc")
    p.add_argument("-o", "--output", required=True, help="New init.nc to create")
    p.add_argument("-t", "--timestep", type=int, default=-1,
                   help="Timestep index (default: -1 = last)")
    p.add_argument("--init-file", required=True, help="Original init.nc")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--no-flip", action="store_true",
                   help="Skip layer reordering")
    p.add_argument("--keep-time", action="store_true",
                   help="Preserve output time value (may crash Init_Spawning "
                        "in current Atlantis builds)")
    p.add_argument("--force", action="store_true",
                   help="Overwrite output if exists")
    args = p.parse_args()
    if not os.path.isfile(args.input):
        p.error(f"Not found: {args.input}")
    if os.path.isfile(args.output) and not args.force:
        p.error(f"Exists: {args.output} (use --force)")
    if not os.path.isfile(args.init_file):
        p.error(f"Not found: {args.init_file}")
    return args


def main():
    args = parse_args()

    # ================================================================
    # Step 1: Copy the original init file
    # ================================================================
    print(f"Copying: {args.init_file} -> {args.output}")
    if os.path.isfile(args.output) and args.force:
        os.remove(args.output)
    shutil.copy2(args.init_file, args.output)

    # ================================================================
    # Step 2: Open files
    # ================================================================
    print(f"Output:  {args.input}")
    ds_out = nc.Dataset(args.input, "r")
    ds_new = nc.Dataset(args.output, "r+")

    wcnz = int(ds_out.getncattr("wcnz"))
    sednz = int(ds_out.getncattr("sednz"))
    nz = wcnz + sednz
    nbox = len(ds_out.dimensions["b"])
    nt_out = len(ds_out.dimensions["t"])

    timestep = args.timestep if args.timestep >= 0 else nt_out - 1
    if timestep >= nt_out:
        print(f"ERROR: timestep {timestep} out of range [0, {nt_out - 1}]")
        sys.exit(1)

    t_val = float(ds_out.variables["t"][timestep])
    t_units = (ds_out.variables["t"].getncattr("units")
               if "units" in ds_out.variables["t"].ncattrs() else "")

    print(f"\nGeometry:  {nbox} boxes x {wcnz} wc + {sednz} sed (z={nz})")
    print(f"Timesteps: {nt_out} in output")
    print(f"Extracting: index {timestep}  ->  t = {t_val} ({t_units})")
    t_human = time_to_human(t_val, t_units)
    if t_human:
        print(f"            {t_human}")

    numlayers = get_numlayers(ds_out, timestep, wcnz)

    # ================================================================
    # Step 3: Update time
    # ================================================================
    if args.keep_time:
        ds_new.variables["t"][0] = t_val
        t_days = t_val / 86400 if "seconds" in t_units.lower() else t_val
        print(f"\nTime:  KEPT at {t_val} ({t_days:.0f} days)")
        print("  WARNING: Current Atlantis builds may crash at Init_Spawning")
        print("  with non-zero start time. Use without --keep-time if this happens.")
    else:
        ds_new.variables["t"][0] = 0.0
        t_days = t_val / 86400 if "seconds" in t_units.lower() else t_val
        print(f"\nTime:  RESET to 0 (extracted from day {t_days:.0f})")
        print("  NOTE: Forcing/harvest files must be offset to match.")

    # ================================================================
    # Step 4: Update tracers (t,b,z) with layer flip
    # ================================================================
    stats = {"tracers": 0, "epi": 0, "phys": 0, "skipped": 0}

    print(f"\nUpdating tracers...")
    for varname in ds_new.variables:
        var_new = ds_new.variables[varname]
        if var_new.dimensions != ("t", "b", "z"):
            continue
        bmtype = (var_new.getncattr("bmtype")
                  if "bmtype" in var_new.ncattrs() else "")

        if bmtype == "tracer" and varname in ds_out.variables:
            data = ds_out.variables[varname][timestep, :, :]
            if hasattr(data, "filled"):
                data = data.filled(0.0)
            if not args.no_flip:
                data = flip_layers_output_to_input(data, numlayers, wcnz)
            var_new[0, :, :] = data
            stats["tracers"] += 1
            if args.verbose:
                print(f"  {varname}")

        elif bmtype == "phys" and varname == "dz":
            if varname in ds_out.variables:
                data = ds_out.variables[varname][timestep, :, :]
                if hasattr(data, "filled"):
                    data = data.filled(0.0)
                if not args.no_flip:
                    data = flip_layers_output_to_input(data, numlayers, wcnz)
                var_new[0, :, :] = data
                stats["phys"] += 1
                if args.verbose:
                    print(f"  {varname}  (phys, flipped)")
        else:
            stats["skipped"] += 1

    # ================================================================
    # Step 5: Update epi (t,b) and physical (t,b) variables
    # ================================================================
    phys_tb_update = {"numlayers", "topk", "sedbiodens", "sedbiodepth",
                      "seddetdepth", "sedirrigenh", "sedoxdepth", "sedturbenh"}

    print(f"Updating epi and physical (t,b) variables...")
    for varname in ds_new.variables:
        var_new = ds_new.variables[varname]
        if var_new.dimensions != ("t", "b"):
            continue
        bmtype = (var_new.getncattr("bmtype")
                  if "bmtype" in var_new.ncattrs() else "")

        if bmtype == "epibenthos" and varname in ds_out.variables:
            data = ds_out.variables[varname][timestep, :]
            if hasattr(data, "filled"):
                data = data.filled(0.0)
            var_new[0, :] = data
            stats["epi"] += 1
            if args.verbose:
                print(f"  {varname}  (epi)")

        elif bmtype == "phys" and varname in phys_tb_update:
            if varname in ds_out.variables:
                data = ds_out.variables[varname][timestep, :]
                if hasattr(data, "filled"):
                    data = data.filled(0.0)
                var_new[0, :] = data
                stats["phys"] += 1
                if args.verbose:
                    print(f"  {varname}  (phys)")

    # ================================================================
    # Step 6: Update history
    # ================================================================
    ds_new.setncattr(
        "history",
        f"Restart from {os.path.basename(args.input)} "
        f"timestep {timestep} (t={t_val})"
    )

    # ================================================================
    # Summary
    # ================================================================
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    print(f"  Tracers updated:  {stats['tracers']}")
    print(f"  Epi updated:      {stats['epi']}")
    print(f"  Physics updated:  {stats['phys']}")
    print(f"  Unchanged:        {stats['skipped']}")

    # Verify
    print(f"\nData check — box 0 ({int(numlayers[0])} layers):")
    for vn in ("dz", "Temp"):
        if vn in ds_out.variables and vn in ds_new.variables:
            o = ds_out.variables[vn][timestep, 0, :]
            n = ds_new.variables[vn][0, 0, :]
            if hasattr(o, "filled"): o = o.filled(0.0)
            if hasattr(n, "filled"): n = n.filled(0.0)
            print(f"  {vn:6s} out: {np.array2string(o, precision=4, suppress_small=True)}")
            print(f"  {vn:6s} new: {np.array2string(n, precision=4, suppress_small=True)}")

    shallow = np.where((numlayers > 0) & (numlayers < wcnz))[0]
    if len(shallow) > 0:
        sb = shallow[0]
        print(f"\nData check — box {sb} ({int(numlayers[sb])} layers, shallow):")
        for vn in ("dz", "Temp"):
            if vn in ds_out.variables and vn in ds_new.variables:
                o = ds_out.variables[vn][timestep, sb, :]
                n = ds_new.variables[vn][0, sb, :]
                if hasattr(o, "filled"): o = o.filled(0.0)
                if hasattr(n, "filled"): n = n.filled(0.0)
                print(f"  {vn:6s} out: {np.array2string(o, precision=4, suppress_small=True)}")
                print(f"  {vn:6s} new: {np.array2string(n, precision=4, suppress_small=True)}")

    ds_new.close()
    ds_out.close()

    fsize = os.path.getsize(args.output) / (1024 * 1024)
    print(f"\nFile: {os.path.abspath(args.output)} ({fsize:.1f} MB)")
    if not args.keep_time:
        print(f"\nRun with:  -i {os.path.basename(args.output)} 0")
    else:
        print(f"\nRun with:  -i {os.path.basename(args.output)} 0")
        print(f"           toutstart = {t_days:.0f} day  (in Run.prm)")


if __name__ == "__main__":
    main()
