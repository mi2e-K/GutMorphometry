# GutMorphometry

Fiji/ImageJ macro and Python helper scripts for semi-automated morphometry of
mouse gut Swiss-roll H&E sections.

The tool measures:

- villus height
- villus width
- villus area
- crypt depth
- muscularis thickness

The Fiji macro handles the interactive workflow and ROI collection. Python
scripts perform the image-processing steps used for villus width and
muscularis thickness QC.

## Repository Layout

```text
fiji/
  plugins/
    Gut_Morphometry.ijm
  scripts/Plugins/AutoRun/
    musc_thickness.py
    tune_muscularis.py
    villus_width.py
docs/
  GutMorphometry_Guide.md
requirements.txt
```

## Requirements

- Windows 10/11
- Fiji/ImageJ
- Python 3.8 or later
- Python packages listed in `requirements.txt`

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

## Installation

Copy the files into a Fiji installation:

```text
fiji/plugins/Gut_Morphometry.ijm
  -> [Fiji.app]/plugins/Gut_Morphometry.ijm

fiji/scripts/Plugins/AutoRun/*.py
  -> [Fiji.app]/scripts/Plugins/AutoRun/
```

Then edit `g_python_cmd` near the top of `Gut_Morphometry.ijm` if Fiji cannot
find your Python executable through `PATH`.

## Workflow Summary

Open a gut section image in Fiji and run:

```text
Plugins > Gut Morphometry [G]
```

The expected image filename format is:

```text
[AnimalID]_[Region]_[Section].tif
```

Region codes:

- `U`: upper small intestine
- `M`: mid small intestine
- `L`: lower small intestine
- `C`: colon

Output files are written next to each image, usually in an image-specific
subfolder:

- `*_results.csv`
- `*_coords.jsonl`
- `*_outer_V*_villus.png`
- `*_outer_MU*_muscularis.png`
- `*_overlay.jpg`

See [docs/GutMorphometry_Guide.md](docs/GutMorphometry_Guide.md) for the full
Japanese user guide and algorithm notes.

## Notes

This repository intentionally excludes raw microscopy images and generated
analysis outputs. Keep large `.tif`, `.png`, `.jpg`, `.csv`, and `.jsonl`
result files outside Git unless they are deliberately curated as examples.

## License
MIT
