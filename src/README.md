# BCG-centered group lensing ESD measurement

This C++17/MPI program measures excess surface density (ESD) around the brightest group galaxies (BCGs). It uses group redshift and halo mass from the group catalogs. MPI rank 0 writes the eleven-column measurement table and renders the figures from the in-memory results.

## Build and run

The build uses an MPI C++17 compiler, a CPython development environment with `python3-config --embed` support, NumPy and Matplotlib. The Python development and runtime environment should be the same. `make` embeds the plotting module in `libesd_plot.a` and links it into `esd_bcg_mpi`.

From the source directory, activate the appropriate Python environment and run:

```bash
make -j4
cp esd.example.ini esd.ini
mpirun -np 3 ./esd_bcg_mpi
```

Set `group-list` and `shear-list` to your data lists before launching. The executable takes no arguments. Compiled defaults are collected in [esd_defaults.hpp](esd_defaults.hpp). Rank 0 reads `esd.ini` from the run directory when present. Set `FQ_GGLENS_INI` to select another file, for example `FQ_GGLENS_INI=sample.ini mpirun -np 3 ./esd_bcg_mpi`. Keys omitted from the INI retain their compiled defaults. Paths in the INI are relative to the run directory, and paths within a catalog list are relative to that list. INI values may be quoted. Use `#` or `;` on a separate line for comments.

On pi2, build in a compute allocation with `module load gcc/12.3.0 openmpi/4.1.6-gcc-12.3.0` and `conda activate pysx` (CPython 3.12, NumPy 2.4.3, Matplotlib 3.10.9). Use `make PYTHON_CONFIG="$CONDA_PREFIX/bin/python3-config"` and run `mpirun -np "$SLURM_NTASKS" ./esd_bcg_mpi` in the batch script. Move run-specific settings from older Slurm commands into an INI for each run. The checked local toolchain was GCC 15.2.0, Open MPI 5.0.10, CPython 3.12.13, NumPy 2.5.0 and Matplotlib 3.11.0.

For the local Pixi installation, set `PIXI_BIN` and `PIXI_PREFIX` to the executable and environment directories, then build and run.

```bash
make -j4 CXX="$PIXI_BIN/mpicxx" PYTHON_CONFIG="$PIXI_PREFIX/bin/python3-config"
"$PIXI_BIN/mpirun" -np 3 ./esd_bcg_mpi
```

## Inputs and BCG matching

`[run] group-list` points to a list file. Its first nonblank line must be `# group_cat galaxy_cat member_cat`. Each block then has four nonblank lines: a unique positive integer block number, the group-catalog path, the galaxy-catalog path, and the member-catalog path, in that order. Blank lines and later `#` comments are ignored. `group_catalogs.list` contains a two-block NGC/SGC example. Catalog file names and locations are unrestricted; paths in a list may be absolute or relative to the list file's directory, and may be enclosed in matching quotes. The program streams each group file and retains lenses passing the half-open redshift and mass cuts.

Each numbered block is an independent group-ID namespace. The same numeric group ID may occur in different blocks. For a block with selected lenses, the program streams its listed member file once, takes column 1 as galaxy ID, column 2 as group ID, and column 3 equal to 1 as the BCG indicator. It then streams the listed galaxy file once and takes columns 1–2 (RA, Dec) from the 1-based row numbered by that galaxy ID. Missing or duplicate BCGs within one block are fatal errors. A `group -> member galaxy ID` file without a BCG flag is not used.

`[run] shear-list` takes one catalog path per nonblank line. Quoted paths and `#` comments are accepted. Each catalog begins with one header line. At startup, rank 0 counts the named columns in the first header and maps the used FQ column names case-insensitively after `EXPO_NUM`. The other readers skip the first line and use the broadcast schema. Use one common layout per run, including the order of columns. Extra unused columns can appear in either part of the schema.

For a legacy catalog without `EXPO_NUM`, rank 0 uses `CCD_NUM` as the boundary of the FQ suffix. The per-exposure bar requires `EXPO_NUM`. See [LEGACY_CATALOG.md](LEGACY_CATALOG.md).

Specify one-based positions in the external prefix under `[external_columns]`. These positions are independent of header spelling.

| Key | DES default | Used for |
| --- | ---: | --- |
| `ra`, `dec` | 5, 6 | High-precision sky position in degrees and adjacent duplicate grouping |
| `zp` | 17 | Source photometric redshift |
| `mag` | 11 | Magnitude limits and star-bar fitting |
| `flags-footprint`, `flags-foreground`, `flags-gold`, `ext-mash` | 1, 2, 3, 4 | Corresponding external flag cuts when enabled |

The adjacent-position grouping uses a tolerance of `1e-7` degrees. Mapped numeric fields with NaN/Inf cause the row to be skipped.

| INI key under `[cuts]` | Default | Selection rule |
| --- | ---: | --- |
| `snrfcut`, `snrlow`, `snrhigh` | 0, 0, 0 | Fourier SNR lower bound; positive lower/upper bounds for `h_flux/sqrt(h_area)` |
| `starcut`, `chi2-thresh` | 16, 0.01 | star statistic ≥ cut and row chi2 ≤ cut |
| `imaxcut`, `jmaxcut` | 64, 64 | maxima strictly below cuts |
| `mag-low`, `mag-high` | 0, 25 | inclusive interval for the configured magnitude |
| `zplow`, `zphigh` | 0, 5 | strictly inside photo-z interval |
| `gf-lim` | 0.0015 | absolute field distortion ≤ cut for both components |
| `psf-chi2-mltp`, `r-half-thresh` | 3, 0 | polynomial PSF chi2 ≤ cut; positive half-size threshold enables size cut |
| `ft-cut`, `fg-cut`, `gold-cut`, `ext-cut` | -1, -10, -10, -4 | nonnegative values enable equality cuts (tolerance `1e-3`) |
| `max-dup` | 5 | reject a whole adjacent-position group above this count |
| `star-bar-mltp`, `expo-bar-mltp` | -3, -3 | global bar activates at ≥0; per-exposure bar activates at >0 and takes precedence |
| `flagcut` | 0 | disabled at 0; a positive value optionally enables the flag × significance upper bound |

The flag cut is commented out in the reference reader, so its default is inactive here. The bar settings default to disabled. When enabled, the program follows the reference size–magnitude histogram and sigma-threshold procedure across MPI ranks before pairing sources. The global bar remains available for legacy input. A positive `expo-bar-mltp` requires a unique `expo_num` field and fails immediately on the legacy header, even when no rows pass cuts. The CCD number cannot identify a physical exposure. Rows with nonpositive `h_area` or PSF size are rejected because the size and SNR expressions are undefined there.

`[run] cache` is optional. A valid cache stores the fully selected BCG lens table and avoids rescanning the large group, member, and galaxy catalogs on later runs with the same cuts and input file identities. The key includes ordered block numbers, canonical input paths, sizes, modification timestamps, and lens selection limits. A stale cache, including the previous fixed-name cache format, is regenerated atomically. Store it in a writable directory.

## Measurement and output

Set `[calibration] fd-corrected = true` if the shear catalog already removed the full two-component FD response from `g1/g2`. The reader restores raw Cartesian moments before either measurement pass. `[calibration] self-calibrate = true` fits both components on the union of sources matched to at least one lens. Otherwise the program uses the configured `m1/m2/c1/c2` values, which default to zero. See [SELF_CALIBRATION.md](SELF_CALIBRATION.md) for the exact moment transformation.

The C++ program calculates comoving distances for flat Lambda-CDM, streams local source catalogs in bounded chunks, searches each chunk around BCGs with the same spherical 3D KD-tree and geometry, applies the lens-source redshift cut, and writes tangential and cross-component PDF_SYM inputs into compact radial scratch files. It uses logarithmic radial bins in comoving `h^-1 Mpc`. BCG positions also determine spherical k-means jackknife regions. PDF_SYM uses the original equal-probability bins, sign-test chi-squared search, and weighted jackknife expression. Candidate histograms are batched into MPI reductions to reduce communication. The NFW reference uses an M_180m profile and the selected concentration model.

The result table contains `R_low R_high R_eff n_pair ESD_meas ESD_sigma ESD_cross ESD_cross_sigma ESD_true ratio residual`. Comment lines record the effective calibration and whether catalog FD was restored. For `ESD_results.dat`, the ESD figure defaults to `ESD_results.png`. Self-calibration also produces `ESD_results_FD.png`. Set `[plots] esd-plot` and `fd-plot` to choose other paths. Rank 0 generates both after measurement using the linked plotting library.

## Bounded memory and scratch storage

- `[run] memory-mb` sets the buffer budget per MPI rank in **MiB**, default **512**, minimum 64. It controls the source chunk, total bin write buffers, current-bin RAM cache, pair scan block and exact-quantile sort block. It is a working-buffer budget, not an OS-enforced RSS limit; lens/file metadata, MPI, allocator overhead and filesystem cache are additional. Internal dimension checks reject oversized duplicate/bin/replica/exposure metadata configurations; at most 256 R bins are supported to bound simultaneously open spool files. At 40 ranks/node the default reserves substantial room beneath the 160 GiB/node allocation.
- `[run] scratch-dir` sets the writable scratch parent. Default: `.bcg_scratch` beside the absolute output path. Each MPI invocation creates a unique run directory and a separate rank directory. A shared filesystem works, and identically named node-local directories can also be used when provisioned on every allocated node.

The flow first streams quality/MAX_DUP-selected source rows with raw shear moments and optional star-bar metadata into binary scratch files. The global stellar bar is fitted once before lens association; exposure bars use bounded histogram batches but keep the original fallback averaged over all successful exposures. Source rows then pass the fitted threshold and enter bounded KD-tree chunks. Original input-file ownership among ranks is preserved. Adjacent-position groups are completed before a source row is emitted, so chunk boundaries cannot split MAX_DUP decisions.

Every pair is written directly into exactly one R-bin file. Its scratch record contains four **double** PDF inputs and the lens index (40 bytes); log(R), lens presence and region counts are accumulated separately. There is no full-pair vector or second radial copy. Source/bar scratch files are removed after pairing. During measurement only one small bin may be cached, up to half the rank budget; larger bins stay on disk and are scanned through fixed blocks. A completed bin is released and its file removed immediately.

Jackknife uses an omitted-region predicate on those same records, with no g/response sample copies. Initial rank-local full and replica quantiles remain exact: bounded sorting and <=32-way external merges derive all omitted-region order statistics from one sorted component stream. This does not use sampled or approximate quantiles. Original rank-averaged boundaries, occupancy refinement (including full-sample normalization for replicas), zero-residual handling, coarse expansion, threshold-20 shrinking, 200-point quadratic fit and weighted jackknife variance are retained. A grid scan reads a block once for all candidate values, then MPI reduces only the small histograms. Chunk order can change floating-point summation order; bitwise equality is not claimed.

**Disk sizing:** pair files require 40 bytes per pair across all ranks. The observed 9,302,755,780-pair sample therefore needs about **372 GB / 347 GiB** for pair payload, plus source/bar spools and temporary exact-sort runs. Sorting uses sixteen-byte magnitude/region records; intermediate merges can temporarily retain both their inputs and output. Provision at least **600 GB free space** for that sample with this default setup, and additional headroom for larger inputs, other jobs or quotas. This is a sizing estimate, not a measured peak. Space requirements are for the whole job, distributed across ranks. Short writes/disk-full errors abort the measurement with the scratch filename.

Normal completion deletes owned scratch files. Caught failures attempt cleanup of the affected rank; SIGKILL, node failure or job cancellation can leave its unique run directory behind. No previous run directory is reused or automatically deleted. Removing an abandoned run's directory is a separate user action. The program recomputes pairs on a new run; it does not resume from an old scratch spool.

## Spherical pairing correction (2026-10-11)

`esd_pairs.cpp` now uses the great-circle separation of the existing source/lens unit vectors for the final radial cut and bin assignment: `R = chi_l * atan2(|source x lens|, source . lens)`. This matches the KD-tree's spherical chord search. The spin-2 and spin-4 rotations use the direction toward the lens in the source's local east/north frame. The Fourier_Quad moment/response signs, beta, selection cuts and measurement algorithms are unchanged.

The geometry regression injects known spin-2/spin-4 moments on six exact spherical rings, including large separations, RA wrap, a small separation and a near-pole case. It runs the production catalog reader, KD tree and pair spool code, then checks pair counts, bins, effective radius, tangential/cross moments and both responses. Before the correction, the z=0.1, Dec=-60 degree, R=90 h^-1 Mpc ring retained 584/720 pairs. After the correction, all six rings retain 720/720 with recovered moments/responses within 1e-8 in real MPI tests with one and three local ranks.

## Initial ESD search interval (2026-10-12)

`measure_bin` uses the per-bin NFW prediction as the tangential initial center and sets the half-width to `w = max(5 * abs(ESD_NFW), 10)` in h M_sun/pc^2. The initial tangential interval is `[ESD_NFW - w, ESD_NFW + w]` for both the full sample and jackknife replicas. The cross component retains its zero center and shares the half-width, giving `[-w, w]`. Subsequent bracket expansion, shrinking and quadratic fitting are retained.

For example, NFW=0.05 gives w=10 and the tangential interval [-9.95, 10.05]; NFW=30 gives w=150 and [-120, 180]. The spherical pairing correction remains in place.
