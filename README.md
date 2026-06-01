# UMass Amherst Biofluids Lab - Mycelial Network Analysis

**Cooper Richman**  
Supervisor: Dr. Sang Hyun Lee  
University of Massachusetts Amherst - Fall 2025 and Spring 2026

---

## Repository Structure

```
C:.
│   .gitignore
├───.vtk & 3D Rendering
├───Data
│   ├───.nd2 Files
│   ├───.tif Files
│   │   ├───Makeshift Data
│   │   ├───Time Stacks
│   │   └───Z Stacks
│   ├───Cleaned .vtk's
│   └───Raw .vtk's
├───Documentation
├───Postprocessing
├───Skeletonization
│   ├───2D
│   │   └───skeleton_plots
│   └───3D
└───Tracking Endpoints & Junctions
```

---

## Directory Descriptions

### `.vtk & 3D Rendering`
Contains the `.vtk Processing.ipynb` notebook for loading, cleaning, and visualizing triangulated surface mesh files in VTK format. Cells cover interactive 3D viewing, batch vertex welding and largest-component extraction, and saving cleaned meshes. The cleaned outputs are written to `Data/Cleaned .vtk's` for use in downstream skeletonization.

---

### `Data`
All raw and intermediate data files used across the project. Subdivided as follows:

#### `Data/.nd2 Files`
Raw microscopy files from the NIS-Elements imaging system. This is not seen here because it is in `.gitignore`. The primary dataset used in this project is:

> `ChannelFITC,TD_Seq0001-001.nd2`

**Note:** `.nd2` files are listed in `.gitignore` and are not tracked by Git because their file sizes are too large for version control. The filenames are documented here for reference. To reproduce results, create the `Data/.nd2 Files` folder, and place the `.nd2` files in this directory before running any notebooks that reference them.

#### `Data/.tif Files/Makeshift Data`
Manually prepared binary TIFF stacks used for early testing and development. Provided in both black-on-white and white-on-black polarity variants.

#### `Data/.tif Files/Time Stacks`
Multi-frame TIFF stacks where **each frame represents a point in time** — i.e., the stack encodes a time-lapse sequence of mycelial growth. These are the primary input format for the Mycelial Path Viewer & Export Tool. Files include cropped, thresholded, and skeletonized variants of the same growth sequence.

#### `Data/.tif Files/Z Stacks`
Multi-frame TIFF stacks where **each frame represents a slice through 3D space** rather than a moment in time. The stack encodes volumetric depth, with each frame being one Z-plane of a three-dimensional structure. These are produced by the 3D skeletonization pipeline and are used for volumetric inspection and smoothing experiments.

#### `Data/Raw .vtk's`
Unprocessed VTK surface mesh files, one per timeframe. These contain noise and disconnected fragments and should be cleaned using `.vtk Processing.ipynb` before use.

#### `Data/Cleaned .vtk's`
VTK mesh files after vertex welding (tolerance = 1.5) and largest-component extraction, produced by `.vtk Processing.ipynb`. These are the input for the skeletonization notebooks.

---

### `Documentation`
PDF documentation for each major component of the project. Refer to these before using any of the scripts or notebooks.

| File | Covers |
|------|--------|
| `Mycelial Path Viewer and Export Tool Software Documentation.pdf` | Full documentation for the main tracking and export tool |
| `.vtk Processing Documentation.pdf` | `.vtk Processing.ipynb` notebook |
| `2D Skeletonization Documentation.pdf` | `2D Skeletonization.ipynb` notebook |
| `3D Skeletonization Documentation.pdf` | `3D Skeletonization.ipynb` notebook |

---

### `Postprocessing`
Contains `Regraphing Growth.py`, the Python source for the companion regraphing script. This script reads the CSV outputs from the Mycelial Path Viewer & Export Tool and generates per-hypha growth plots with summary statistics per frame.

---

### `Skeletonization`

#### `Skeletonization/2D`
Contains `2D Skeletonization.ipynb`, a notebook with independent cells for generating 2D skeletons from `.vtk` files or `.nd2` files. Also contains `temp_mesh.stl`, a temporary intermediate file written during voxelization (overwritten on each run).

#### `Skeletonization/2D/skeleton_plots`
PNG output images from the batch processing cell in `2D Skeletonization.ipynb,` one labeled skeleton plot per frame (frames 0–15).

#### `Skeletonization/3D`
Contains `3D Skeletonization.ipynb`, a notebook for computing true volumetric skeletons using the Lee 3D thinning algorithm, as well as hybrid 2D-to-3D mapping and smoothing experiments.

---

### `Tracking Endpoints & Junctions`
Contains `Mycelial Path Viewer Export Tool.py`, the full Python source code for the main GUI application. See `Documentation/Mycelial Path Viewer and Export Tool Software Documentation.pdf` for complete usage instructions.
