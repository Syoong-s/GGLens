# Developer notes

The executable measures one group sample around one selected center and can fit that measurement immediately or in a later run. C++ owns distributed catalog ingest, source selection, pair creation, ESD estimation, and covariance. The Python modules handle the fit and its output data.

| Reference | Scope |
| --- | --- |
| [Catalogs and calibration](CATALOGS_AND_CALIBRATION.md) | DR9 joins, centroid definitions, shear-column resolution, source selection, and field distortion |
| [ESD and covariance](MEASUREMENT.md) | Pair geometry, response conventions, PDF_SYM, sky regions, and radial covariance |
| [Model, sampling, and products](FITTING_AND_OUTPUTS.md) | Halo terms, power spectrum, likelihood, posterior draws, and figure generation |
| [Publishing releases](RELEASING.md) | Source archive contents and version tags |

## Execution order

1. Rank 0 reads the INI selected by the optional positional argument, or looks for `esd.ini` in the launch directory. It broadcasts the INI text. All ranks apply it to the compiled defaults and validate the settings for the selected stage. A fit-only run reads and broadcasts the saved fit input and resumes at step 5.
2. Rank 0 reads the shear list and resolves the source fields from the first header. Resolve lens column positions from the first catalog of each type, then load the selected centers or restore a matching lens cache. Every rank receives the same ordered lens table.
3. Source workers apply catalog and calibration selections. They write pair records by radial bin. The participating lenses set the sky regions and the fixed fit inputs.
4. Measure full-sample and deleted-region ESD profiles. Rank 0 writes these results and their covariance, then publishes the versioned fit input. A measurement-only run finishes here.
5. Rank 0 computes the matter correlation and broadcasts it. Assigned ranks run independent ensembles. Rank 0 gathers the chains and publishes the fit data and figure.

## Code locations

| Stage | Implementation |
| --- | --- |
| Program entry and configuration | [esd_bcg_mpi.cpp](../code/src/esd_bcg_mpi.cpp), [esd_config.cpp](../code/src/esd_config.cpp), [esd_defaults.hpp](../code/include/esd_defaults.hpp) |
| Lens selection and caches | [esd_lenses.cpp](../code/src/esd_lenses.cpp), [bcg_dr9_manifest.cpp](../code/src/bcg_dr9_manifest.cpp), [bcg_dr9.cpp](../code/src/bcg_dr9.cpp) |
| Weighted centers | [bcg_dr9_centroids.cpp](../code/src/bcg_dr9_centroids.cpp), [bcg_dr9_join.cpp](../code/src/bcg_dr9_join.cpp), [bcg_dr9_rows.hpp](../code/include/bcg_dr9_rows.hpp) |
| Source selection and calibration | [bcg_shear_catalog.cpp](../code/src/bcg_shear_catalog.cpp), [bcg_shear_cuts.cpp](../code/src/bcg_shear_cuts.cpp), [bcg_star_bar.cpp](../code/src/bcg_star_bar.cpp), [selfcalibration.cpp](../code/src/selfcalibration.cpp) |
| Pair storage and measurement | [esd_pairs.cpp](../code/src/esd_pairs.cpp), [bcg_pair_store.cpp](../code/src/bcg_pair_store.cpp), [esd_measure.cpp](../code/src/esd_measure.cpp) |
| Sky regions and covariance | [esd_regions.cpp](../code/src/esd_regions.cpp), [esd_covariance.cpp](../code/src/esd_covariance.cpp) |
| Saved fit input | [esd_fit_input.cpp](../code/src/esd_fit_input.cpp) |
| Fit and plotting bridge | [esd_fit_mpi.cpp](../code/src/esd_fit_mpi.cpp), [esd_fit_backend.cpp](../code/src/esd_fit_backend.cpp), [esd_fit.py](../code/script/esd_fit.py) |

## Build and runtime

The Makefile in `code/` compiles C++17 with `mpicxx`. It uses `PYTHON_CONFIG --includes` and `PYTHON_CONFIG --ldflags --embed` from the same Python environment. The build converts `esd_nfw.py`, `esd_two_halo.py`, `esd_archive.py`, `esd_profiles.py`, `esd_plot.py`, and `esd_fit.py` into generated includes and links them through `libesd_fit.a`. Edit the Python files under `code/script/` and rebuild the executable to update the embedded modules.

The fit bridge initializes isolated CPython only on ranks assigned an ensemble. It disables Python bytecode writes and signal-handler installation. The bridge adds `[fit] python-path` to the embedded search path when set. Use that setting for dependencies outside the linked environment.

The embedded Python interpreter and its compiled extensions share the process's C++ runtime. On Linux, the default `CXX_RUNTIME=auto` looks for `libstdc++.so.6` under `lib/` or `lib64/` of the prefix returned by `PYTHON_CONFIG --exec-prefix`. When `patchelf` and `ldd` are available, the build tries an absolute dependency on that one library in a temporary copy of the executable. It publishes the patched copy after checking the linked dependencies and verifying that MPI and `libpython` resolve to the same paths. If automatic selection is unavailable or the loader rejects the candidate, the build reports the reason and publishes the original link result.

The default `PATCHELF=auto` searches `PATH`, the selected Python environment's `bin/` directory, and Conda base. If the selected environment is under a Conda installation's `envs/`, the build checks that installation directly. For other prefixes, it can query `conda info --base`. Set `PATCHELF` to a specific executable to override discovery.

Use `CXX_RUNTIME=system` to use normal linker resolution or set `CXX_RUNTIME` to an absolute `libstdc++.so.6` path to require a particular runtime. A build with a pinned runtime uses the selected library path at startup. Rebuild after moving the Python environment to another location. Use `test-stages` to check imports, including dependencies loaded only during fitting.

The cluster smoke tests built with GCC 12.3.0 and Open MPI 4.1.6. The runtime was Python 3.12.13, NumPy 2.4.3, SciPy 1.17.1, Matplotlib 3.10.9, CAMB 1.6.5, and emcee 3.1.6. `cumulative_simpson` requires SciPy 1.12 or later.

```sh
make -C code
mpirun -np "$SLURM_NTASKS" ./code/esd_bcg_mpi "/path/to/esd.ini"
```

For builds using normal linker resolution, a mixed environment that reports `CXXABI_*` or `GLIBCXX_*` version errors can pass its compatible runtime to MPI ranks with `-x LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6"`. Choose `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `MKL_NUM_THREADS=1` for MPI tasks allocated one CPU each.

## Validation cases

Run `make -C code test-catalogs` with an MPI installation that can launch four ranks. Set `MPIEXEC` and `MPIEXEC_FLAGS` for the local launcher. The test compares BCG, LWC, and NWC centers from DR9, reordered CSV, and catalogs that identify objects by data row number. It also exercises minimal NWC inputs, startup column errors, and cache invalidation.

Run `make -C code test-stages` with the fit dependencies and two MPI ranks. It checks the lossless archive round trip, combined and split runs, offline regeneration of the displays, and older fit inputs. Set `FQ_GGLENS_TEST_PYTHON_PATH` if the embedded interpreter needs an additional package directory.

| Area | Reference check |
| --- | --- |
| Catalog matching | Check the same Galaxy path in multiple blocks, overlapping Group IDs across blocks, Galaxy ID zero, sparse 64-bit IDs, and rows crossing byte and batch boundaries. Exercise missing and duplicate matches and invalid LWC luminosities. |
| Geometry | Use known member centroids on either side of RA zero and known source rings. Verify exact pair counts and effective radii, and compare the BCG measurement with an established baseline. |
| Measurement | Check optimized sign histograms against direct trial counts, including edge ties and omitted regions. Recompute the full covariance from the saved replicas. |
| Fit | Compare the halo projections to direct quadrature. Check that the minimum retained chain chi-squared agrees with the best saved model evaluated against the full covariance. |
| Products | Recompute the curve percentiles from the indexed draws. Verify that the plot can be regenerated pixel for pixel in the same environment and that a mismatched fit ID is rejected. |
