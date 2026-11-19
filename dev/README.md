# Developer notes

The executable measures one group sample around one selected center and then fits that measurement in the same run. C++ owns distributed catalog ingest, source selection, pair creation, ESD estimation, and covariance. The Python modules handle the fit and its output data.

| Reference | Scope |
| --- | --- |
| [Catalogs and calibration](CATALOGS_AND_CALIBRATION.md) | DR9 joins, centroid definitions, shear-column resolution, source selection, and field distortion |
| [ESD and covariance](MEASUREMENT.md) | Pair geometry, response conventions, PDF_SYM, sky regions, and radial covariance |
| [Model, sampling, and products](FITTING_AND_OUTPUTS.md) | Halo terms, power spectrum, likelihood, posterior draws, and figure generation |
| [Publishing releases](RELEASING.md) | Source archive contents and version tags |

## Execution order

1. Rank 0 reads the INI selected by the optional positional argument, or looks for `esd.ini` in the launch directory. It broadcasts the INI text. All ranks apply it to the compiled defaults and validate the result. Rank 0 reads the shear list and resolves the source fields from the first header.
2. Resolve the selected DR9 lens centers or restore a matching lens cache. Every rank receives the same ordered lens table.
3. Source workers apply catalog and calibration selections. They write pair records by radial bin. The participating lenses set the sky regions and the fixed fit inputs.
4. Measure full-sample and deleted-region ESD profiles. Rank 0 writes these results and their covariance before sampling starts.
5. Rank 0 computes the matter correlation and broadcasts it. Assigned ranks run independent ensembles. Rank 0 gathers the chains and publishes the fit data and figure.

## Code locations

| Stage | Implementation |
| --- | --- |
| Program entry and configuration | [esd_bcg_mpi.cpp](../src/esd_bcg_mpi.cpp), [esd_config.cpp](../src/esd_config.cpp), [esd_defaults.hpp](../include/esd_defaults.hpp) |
| Lens selection and caches | [esd_lenses.cpp](../src/esd_lenses.cpp), [bcg_dr9_manifest.cpp](../src/bcg_dr9_manifest.cpp), [bcg_dr9.cpp](../src/bcg_dr9.cpp) |
| Weighted centers | [bcg_dr9_centroids.cpp](../src/bcg_dr9_centroids.cpp), [bcg_dr9_join.cpp](../src/bcg_dr9_join.cpp), [bcg_dr9_rows.hpp](../include/bcg_dr9_rows.hpp) |
| Source selection and calibration | [bcg_shear_catalog.cpp](../src/bcg_shear_catalog.cpp), [bcg_shear_cuts.cpp](../src/bcg_shear_cuts.cpp), [bcg_star_bar.cpp](../src/bcg_star_bar.cpp), [selfcalibration.cpp](../src/selfcalibration.cpp) |
| Pair storage and measurement | [esd_pairs.cpp](../src/esd_pairs.cpp), [bcg_pair_store.cpp](../src/bcg_pair_store.cpp), [esd_measure.cpp](../src/esd_measure.cpp) |
| Sky regions and covariance | [esd_regions.cpp](../src/esd_regions.cpp), [esd_covariance.cpp](../src/esd_covariance.cpp) |
| Fit and plotting bridge | [esd_fit_mpi.cpp](../src/esd_fit_mpi.cpp), [esd_fit_backend.cpp](../src/esd_fit_backend.cpp), [esd_fit.py](../script/esd_fit.py) |

## Build and runtime

The Makefile compiles C++17 with `mpicxx`. It uses `PYTHON_CONFIG --includes` and `PYTHON_CONFIG --ldflags --embed` from the same Python environment. The build converts `script/esd_nfw.py`, `esd_two_halo.py`, `esd_profiles.py`, `esd_plot.py`, and `esd_fit.py` into generated includes and links them through `libesd_fit.a`. Edit the Python files under `script/` and rebuild the executable to update the embedded modules.

The fit bridge initializes isolated CPython only on ranks assigned an ensemble. It disables Python bytecode writes and signal-handler installation. The bridge adds `[fit] python-path` to the embedded search path when set. Use that setting for dependencies outside the linked environment.

The cluster smoke tests built with GCC 12.3.0 and Open MPI 4.1.6. The runtime was Python 3.12.13, NumPy 2.4.3, SciPy 1.17.1, Matplotlib 3.10.9, CAMB 1.6.5, and emcee 3.1.6. `cumulative_simpson` requires SciPy 1.12 or later. In the tested mixed compiler and Conda environment, workers preloaded the Conda C++ runtime to satisfy the SciPy extensions' ABI.

```sh
make -j 4 PYTHON_CONFIG="$(command -v python3-config)" CPPFLAGS=-DOMPI_SKIP_MPICXX
mpirun -np "$SLURM_NTASKS" -x LD_PRELOAD="$CONDA_PREFIX/lib/libstdc++.so.6" ./esd_bcg_mpi "/path/to/esd.ini"
```

Apply the preload to a mixed environment that needs it. Choose `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and `MKL_NUM_THREADS=1` for MPI tasks allocated one CPU each.

## Validation cases

| Area | Reference check |
| --- | --- |
| Catalog matching | Check the same Galaxy path in multiple blocks, overlapping Group IDs across blocks, Galaxy ID zero, sparse 64-bit IDs, and rows crossing byte and batch boundaries. Exercise missing and duplicate matches and invalid LWC luminosities. |
| Geometry | Use known member centroids on either side of RA zero and known source rings. Verify exact pair counts and effective radii, and compare the BCG measurement with an established baseline. |
| Measurement | Check optimized sign histograms against direct trial counts, including edge ties and omitted regions. Recompute the full covariance from the saved replicas. |
| Fit | Compare the halo projections to direct quadrature. Check that the minimum retained chain chi-squared agrees with the best saved model evaluated against the full covariance. |
| Products | Recompute the curve percentiles from the indexed draws. Verify that the plot can be regenerated pixel for pixel in the same environment and that a mismatched fit ID is rejected. |
