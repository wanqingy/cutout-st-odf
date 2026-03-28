# cutout-st-odf

A self-contained notebook pipeline for 3-D fluorescence microscopy that:

1. downloads a small volume cutout from a Neuroglancer link (or loads a local `.tif`),
2. computes voxel-wise orientation vectors with a 3-D structure tensor, and
3. renders a direction-encoded ODF (orientation distribution function) glyph.

```
Step 1 – Cutout          →   Step 2 – Structure Tensor   →   Step 3 – ODF Glyph
 zarr store / local tif       st_odf.StructureTensor           odf_from_vectors_simple
```

## Prerequisites

| Dependency | Notes |
|---|---|
| Python ≥ 3.12 | |
| `numpy`, `scipy` | Numerics |
| `pyvista` | Off-screen glyph rendering |
| `matplotlib`, `Pillow` | Inline display |
| `tifffile` | Reading / writing `.tif` cutouts |
| `ac_data_tools` *(optional)* | Only needed for Option A (cutout from Neuroglancer) |

## Install

```bash
git clone https://github.com/wanqingy/cutout-st-odf.git
cd cutout-st-odf
pip install -e .
```

`ac_data_tools` is an internal Allen Institute package. If you have access, install it
separately; otherwise use **Option B** in Step 1 to load a local `.tif` file:

```bash
pip install git+https://github.com/AllenInstitute/ac_data_tools.git
```

## Quick start

1. Add your dataset entry to `datasets_config.json` (see below).
2. Open `notebooks/cutout_st_odf_glyph.ipynb`.
3. Set `DATASET_ID` and (for Option A) paste a Neuroglancer share link into `NG_LINK`.
4. Run all cells.

## Dataset config (`datasets_config.json`)

Each entry in the `"datasets"` array describes one imaging dataset:

```json
{
  "datasets": [
    {
      "id": "MY_DATASET",
      "name": "My Dataset Name",
      "anatomical_axes": ["A-P", "D-V", "M-L"],
      "anatomical_rotation": [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
      ],
      "dec_rgb_order": [0, 1, 2]
    }
  ]
}
```

| Field | Description |
|---|---|
| `id` | Short identifier used in the notebook (`DATASET_ID`) |
| `anatomical_axes` | Human-readable labels for the three output axes after rotation |
| `anatomical_rotation` | 3×3 matrix **R** such that `v_anat = R @ v_image` |
| `dec_rgb_order` | `[r_axis, g_axis, b_axis]` — indices into the rotated vector for DEC colouring (convention: horizontal=R, vertical=G, depth=B) |

Use the 3×3 identity matrix for `anatomical_rotation` if no rotation is needed.

## Axis conventions

| Stage | Axis order | Notes |
|---|---|---|
| zarr store / `StructureTensor` input | **(z, y, x)** | Standard image convention |
| `cutout_from_NG` return value | **(x, y, z)** | `ac_data_tools` transposes internally — assign directly to `img` |
| Structure-tensor eigenvectors | **(v_z, v_y, v_x)** | Matches image axis order |
| Anatomical space | **R @ v_image** | Applied via `anatomical_rotation` from config |

## ODF notes

- **Antipodal symmetry** — Structure-tensor eigenvectors have arbitrary sign.
  The pipeline concatenates `[v, -v]` before calling `odf_from_vectors_simple`
  so that both hemispheres of the sphere are populated.
- **Subsampling** — Large volumes are randomly subsampled to `MAX_VECTORS`
  (default 125 000) before ODF computation to keep runtime short.
- **GFA** — The Generalized Fractional Anisotropy scalar is printed after ODF
  computation; values close to 1 indicate a strongly peaked distribution.

## Package API (`src/st_odf/`)

| Symbol | Description |
|---|---|
| `StructureTensor` | Compute 3-D structure tensor and principal eigenvectors for a `(Z, Y, X)` volume |
| `make_sphere` | Fibonacci-sampled unit sphere mesh (vertices + ConvexHull faces) |
| `odf_from_vectors_simple` | Build a smooth ODF on a sphere from orientation vectors (`vonmises`, `kde`, or `hist` kernel) |
| `generalized_fractional_anisotropy` | GFA scalar from an ODF array |
| `apply_anatomical_rotation` | Rotate sphere vertices by a 3×3 matrix (`R @ v`) |

## License

TBD
