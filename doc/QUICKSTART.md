# Quick start

## 1. Prepare the input lists

Prepare a numbered list of DR9 Group, Galaxy, and Member catalogs and a separate list of shear catalogs. Use the same column order in every shear file for the run. [Input formats](INPUTS.md)

## 2. Create the run configuration

Save a configuration file such as `esd.ini`. The example below uses a shear header that starts with `ra dec zphot mag` and ends with exposure quality. Set the column positions to match your files. [Shear header example](INPUTS.md#shear-catalogs)

```ini
[run]
group-list = group_catalogs.list
shear-list = shear_catalogs.list
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

Omitted settings use their defaults. The default group sample has redshifts from `0.3` up to but excluding `0.4` and catalog log masses from `13.0` up to but excluding `14.2`. Set `z-min`, `z-max`, `logm-min`, and `logm-max` in `[measurement]` to choose another sample. Set `[calibration] fd-corrected = true` if field distortion was already removed from the input shear values.

The default measurement uses 12 radial bins from `0.1` to `100.0` comoving Mpc divided by `h`, with automatic calibration enabled. The values set in the full INI example are the program defaults. To turn off a default cut, set its value below `-999.0`. [Full INI example](../esd.example.ini) · [Configuration guide](CONFIGURATION.md)

## 3. Launch the program

If you need to build the program, activate the C++17, MPI, and Python environments provided for your system. The Python environment needs NumPy, SciPy 1.12 or later, Matplotlib, CAMB, and emcee. From the repository root, run

```sh
make -j 4 PYTHON_CONFIG="$(command -v python3-config)"
```

Launch from the run directory and pass the path to your INI file. Use your allocated task count and the path to your build. On Slurm, use `$SLURM_NTASKS` for the task count.

```sh
mpirun -np "$SLURM_NTASKS" /path/to/FQ_GGLens/esd_bcg_mpi "/path/to/esd.ini"
```

Omit the INI argument to use `esd.ini` in the launch directory. If that file is absent, the program uses its built-in defaults. Relative paths inside the INI start from the launch directory.

The figure and companion files are written beside the path set by `output`. Use a distinct output name for each sample or center you want to keep. [Result files](RESULTS.md)
