# FQ GGLens

FQ GGLens measures the excess surface density (ESD) of galaxy groups from configurable lens catalogs and Fourier Quad shear catalogs. The default lens layout matches DESI Legacy Imaging Surveys DR9. Choose the brightest central galaxy (BCG), the luminosity-weighted center (LWC), or the number-weighted center (NWC) for each run.

A run saves the measured ESD and its uncertainties, the best-fit profile, the range of fitted profiles, and a figure. The tables contain the values needed to plot the results again.

| Guide | What it covers |
| --- | --- |
| [Quick start](doc/QUICKSTART.md) | Build, configure, and launch a first run |
| [Input catalogs](doc/INPUTS.md) | Lens catalog layouts and the required shear columns |
| [Configuration](doc/CONFIGURATION.md) | Center choice, sample limits, and optional settings |
| [Results](doc/RESULTS.md) | Output filenames, units, and table columns |

Download `FQ_GGLens-v<version>.tar.gz` from GitHub Releases. Run `make -C code` from the extracted archive or repository root to build `code/esd_bcg_mpi`. Open the matching tag for these guides and the example INI. [Example configuration](esd.example.ini)

FQ GGLens is distributed under the [MIT License](LICENSE).
