# cutout-st-odf

A self-contained notebook pipeline that downloads a small volume cutout from a
Neuroglancer link, computes structure-tensor orientation vectors,
and renders a orientation-encoded ODF (orientation distribution function) glyph.

```
Step 1 – Cutout     →   Step 2 – Structure Tensor   →   Step 3 – ODF Glyph
 zarr / local tif        fiberorient.StructureTensor      odf_from_vectors_simple
```

## Prerequisites

| Dependency | Notes |
|---|---|
| Python ≥ 3.10 | Tested on 3.12 |
| `fiberorient` | Structure-tensor fitting |
| `pyvista` | Off-screen glyph rendering |
| `matplotlib` | Inline display |
| `numpy`, `scipy` | Numerics |
| `ac_data_tools` *(optional)* | Only needed for Option A (cutout from Neuroglancer) |

## Install

```bash
git clone https://github.com/wanqingy/cutout-st-odf.git
cd cutout-st-odf
pip install -e .
```

`ac_data_tools` is a private/internal package. If you have access, install it
separately; otherwise use **Option B** in Step 1 to load a local `.tif` file.

## Quick start

1. Copy `datasets_config.json` and fill in your dataset entry.
2. Open `notebooks/cutout_st_odf_glyph.ipynb`.
3. Set `DATASET_ID` and (for Option A) paste a Neuroglancer share link into
   `NG_LINK`.
4. Run all cells.

## Dataset config

```json
{
  "MY_DATASET": {
    "anatomical_rotation": [[1,0,0],[0,1,0],[0,0,1]]
  }
}
```

`anatomical_rotation` is a 3×3 rotation matrix that maps structure-tensor
eigenvectors from image space into anatomical (e.g. LPS) space before the ODF
is computed. Use the 3×3 identity matrix if no rotation is needed.

## Axis conventions

| Layer | Axis order | Notes |
|---|---|---|
| zarr store | (z, y, x) | Neuroglancer convention |
| `cutout_from_NG` return value | **(x, y, z)** | ac_data_tools transposes internally |
| `fiberorient.StructureTensor` input | **(z, y, x)** | Fix: `img = cutout.transpose(2,1,0)` |
| Middle-plane display (grayscale) | transpose with `.T` | matches Neuroglancer view |
| Middle-plane display (RGB) | `.transpose(1,0,2)` | spatial-only; preserves colour axis |

## ODF notes

- **Antipodal symmetry** — Structure-tensor eigenvectors have arbitrary sign.
  The pipeline doubles the vector set (`np.concatenate([v, -v])`) before
  calling `odf_from_vectors_simple` so that both hemispheres are populated.
- **Subsampling** — Large volumes are subsampled to `MAX_VECTORS = 125 000`
  before ODF computation to keep runtime short.

## Extracted utilities (`src/st_odf/utils.py`)

| Function | Description |
|---|---|
| `odf_from_vectors_simple` | Builds an ODF on a sphere from a set of unit vectors |
| `generalized_fractional_anisotropy` | GFA scalar from an ODF array |
| `apply_anatomical_rotation` | Rotates sphere vertices by a 3×3 matrix |

## License

TBD
