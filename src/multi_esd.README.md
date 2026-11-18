# Multiple ESD tangential ratios

`multi_esd.py` compares the tangential ESD from multiple completed runs. Its ratio-axis and horizontal-offset helpers are shared with the embedded plotter in `esd_plot.py`. Keep the two files together when running this offline comparison.

## Interpreter and dependencies

- Local check: Python 3.12.13, NumPy 2.5.0, Matplotlib 3.11.0.
- Cluster target: `conda activate pysx`; Python 3.12.13, NumPy 2.4.3,
  Matplotlib 3.10.9 (existing environment, checked using package metadata).
- No compilation or installation is required. Uses Python standard library,
  NumPy and Matplotlib; the noninteractive Agg backend is selected explicitly.

## Run from the source directory

Activate `pysx` on the cluster, then use this portable Linux command:

```bash
python multi_esd.py --input result_A/ESD_results.dat --label "Sample A" --input result_B/ESD_results.dat --label "Sample B" --output multi_esd.png
```

Equivalently, supply the ordered lists together:

```bash
python multi_esd.py --input result_A/ESD_results.dat result_B/ESD_results.dat --label "Sample A" "Sample B" --output multi_esd.png
```
```bash
python multi_esd.py --input /lustre/home/acct-phyzj/share/DES/testy/cluster_lensing/mode_bcg/result-star-bar/ESD_results.dat /lustre/home/acct-phyzj/share/DES/testy/cluster_lensing/mode_bcg/result-star-bar7/ESD_results.dat /lustre/home/acct-phyzj/share/DES/testy/cluster_lensing/mode_bcg/result-lt20/ESD_results.dat --label "2.5-Sigma" "7-Sigma" "LT20" --output multi_esd.png
```

python multi_esd.py --input /lustre/home/acct-phyzj/share/DES/testy/cluster_lensing/mode_bcgo/yband/snr15/ESD_results.dat /lustre/home/acct-phyzj/share/DES/testy/cluster_lensing/mode_bcgo/yband/snr15/ESD_results_m.dat /lustre/home/acct-phyzj/share/DES/testy/cluster_lensing/mode_bcgo/yband/snr40/ESD_results.dat /lustre/home/acct-phyzj/share/DES/testy/cluster_lensing/mode_bcgo/yband/snr40/ESD_results_m.dat --label "SNR15" "SNR15-Calib" "SNR40" "SNR40-Calib" --output multi_esd_yband.png --x-spread 0.08

Every input needs a nonempty label; input and label counts must match. Quote
labels containing spaces. Inputs, outputs, labels and the number of inputs
are supplied by the caller; no dataset paths are hardcoded. The output parent
directory must already exist. As with the original script, saving overwrites
an existing output image. `--title` overrides the figure title.

## Data and scale

Column numbering is one-based: radius is column 3, measured ESD column 5,
ESD sigma column 6, NFW theory column 9 and tangential ratio column 10.
The original ratio is preserved; the uncertainty is sigma / abs(theory).
Radius must be positive, and radius, measured ESD, theory and ratio must be
finite; zero theory and invalid rows are excluded with warnings. Invalid
uncertainties keep the ratio point without an error bar and produce a warning.
An empty/malformed table or an input with no valid ratio points is an error.
Each series is sorted by radius for display; no bin alignment is required.

The horizontal axis is logarithmic. The vertical display coordinate is
`asinh((ratio - 1) / width)`, with a default width of 0.1. This is approximately
linear near 1 and logarithmically compressed far from 1, including negative
ratios. Tick labels always show the original ratio values. Ratios and error
bar endpoints are transformed together; error bars can look asymmetric.
The dashed reference line is ratio=1.

- `--linear-width 0.05`: make the central finely resolved region narrower.
- `--ylim 0.5 1.5`: optional limits in original ratio units. Without this
  option, limits include all series, error bars and the ratio=1 reference line.
- Legend labels, colors and marker shapes distinguish the inputs.
- `--x-spread 0.06` (default): slightly separate dataset points and their error
  bars horizontally. The full spread is 0.06 decades in log10(R), centered
  on the original radius; the outermost series use factors approximately
  0.933 and 1.072. For multiple inputs, offsets are evenly spaced across this
  fixed span, so adding inputs does not expand it. A single input has no
  offset. `--x-spread 0` disables offsets. This changes display coordinates
  only, never the radius column, measured ratios, uncertainties or input files.

## Automatic single-run figure

The MPI executable creates the ESD and tangential-ratio panels at the end of each measurement. Set `[plots] x-spread`, `linear-width`, `ratio-y-min` and `ratio-y-max` in the measurement INI to change its scale or display offsets. See [README.md](README.md).

## Validation performed

The comparison tool has been checked locally with synthetic result tables for the two-input ratio, scale inverse and display offsets. The shared helpers are also exercised by the automatic figures in the MPI smoke runs.
