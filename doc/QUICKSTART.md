# Quick start

## 1. Prepare the input lists

Prepare a numbered list of Group, Galaxy, and Member catalogs and a separate list of shear catalogs. Use the same column order in all files of a given type. The lens column defaults match DR9. [Input formats](INPUTS.md)

## 2. Create the run configuration

Save a configuration file such as `esd.ini`. The example below uses a shear header that starts with `ra dec zphot mag` and ends with exposure quality. Set the column positions to match your files. [Shear header example](INPUTS.md#shear-catalogs)

```ini
[run]
group-list = group_catalogs.list
shear-list = shear_catalogs.list
output-dir = results
output = ESD_results.dat

[measurement]
center = bcg

[external_columns]
ra = 1
dec = 2
zp = 3
mag = 4
expo-chi2 = -1
```

Omitted settings use their defaults. For a different lens layout, set positions in `[group_columns]`, `[galaxy_columns]`, and `[member_columns]`. The default group sample has redshifts from `0.3` up to but excluding `0.4` and catalog log masses from `13.0` up to but excluding `14.2`. Set `z-min`, `z-max`, `logm-min`, and `logm-max` in `[measurement]` to choose another sample. Set `[calibration] fd-corrected = true` if field distortion was already removed from the input shear values. [Lens column settings](CONFIGURATION.md#group-galaxy-and-member-columns)

The default measurement uses 12 radial bins from `0.1` to `100.0` comoving Mpc divided by `h`, with automatic calibration enabled. To turn off a default cut, set its value below `-999.0`. [Full INI example](../code/esd.example.ini) · [Configuration guide](CONFIGURATION.md)

## 3. Launch the program

If you need to build the program, activate the C++17, MPI, and Python environments provided for your system. The Python environment needs NumPy, SciPy 1.12 or later, Matplotlib, CAMB, and emcee. Build from the repository or extracted archive root.

```sh
make -C code
```

The build uses `mpicxx` and `python3-config` from `PATH`. On Linux, it can pin `libstdc++.so.6` from the selected Python environment. It looks for `patchelf` in `PATH`, that environment, and Conda base. [C++ runtime selection](../dev/README.md#build-and-runtime)

Launch from the run directory and pass the path to your INI file. Use your allocated task count and the path to your build. On Slurm, use `$SLURM_NTASKS` for the task count.

```sh
mpirun -np "$SLURM_NTASKS" /path/to/FQ_GGLens/code/esd_bcg_mpi "/path/to/esd.ini"
```

Omit the INI argument to use `esd.ini` in the launch directory. If that file is absent, the program uses its built-in defaults. Relative paths inside the INI start from the launch directory.

With the example settings, the stage archives go in `results/measurement/` and `results/fit/`. The final ESD figure, ESD table and posterior figure go in `results/`. Use a distinct `output` name or `output-dir` for each sample or fit you want to keep. [Result files](RESULTS.md)

## 4. Run the stages separately

Use `--stage` to override `[run] stage` for a launch. The fit command reads the measurement snapshot. The pull command reads both archives and can be rerun when the catalogs are offline. It needs NumPy and Matplotlib.

```sh
mpirun -np "$SLURM_NTASKS" /path/to/FQ_GGLens/code/esd_bcg_mpi --stage measure "/path/to/esd.ini"
mpirun -np "$SLURM_NTASKS" /path/to/FQ_GGLens/code/esd_bcg_mpi --stage fit "/path/to/esd.ini"
mpirun -np 1 /path/to/FQ_GGLens/code/esd_bcg_mpi --stage pull "/path/to/esd.ini"
```

To reuse the measurement in another output directory, make a fit-only INI and launch it with the same executable.

```ini
[run]
stage = fit
output-dir = refits/bcg

[fit]
input = results/measurement/ESD_results_measurement.bin
```

The saved snapshot supplies the ESD, jackknife estimates, center, lens statistics, and cosmology. Set sampling changes in `[fit]` and figure changes in `[plots]`. The new fit records the snapshot location so `pull` can use it from the refit directory. [Fit input settings](CONFIGURATION.md#fit-and-figure)
