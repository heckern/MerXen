# Section alignment *(planned)*

> **Status: not yet implemented.** The stub in
> [src/merxen/alignment/register.py:26](../../src/merxen/alignment/register.py#L26)
> raises `NotImplementedError` and the Nextflow pipeline does **not** invoke
> this stage today.

## Intent

Adjacent tissue sections are run on different platforms, so even after
transforming both to microns they do not share a coordinate system: the
slices are physically offset, and small elastic deformations accumulate during
sample preparation. Comparative analysis at single-cell resolution needs both
platforms co-registered into a shared reference space.

The planned stage will fit between QC and COMPARE in the workflow graph
(see [Pipeline architecture](../pipeline.md#stage-graph)) and produce
transforms that map both sections into a common frame.

## Planned interface

Two public symbols already exist as placeholders in
[src/merxen/alignment/register.py](../../src/merxen/alignment/register.py):

```python
@dataclass(frozen=True)
class TransformResult:
    merscope_to_common: Any      # transform object
    xenium_to_common: Any        # transform object
    metadata: dict[str, Any]     # implementation-specific information


def register_pair(
    merscope_sdata: Any,
    xenium_sdata: Any,
    config: Any,
) -> TransformResult:
    ...
```

Future Nextflow module (name placeholder): `ALIGN`.
CLI entry-point placeholder: `merxen align --config align_config.json`.

## What will likely change in the pipeline

- A new Nextflow module `workflows/modules/alignment.nf` invoked after
  per-dataset QC, before `COMPARE`.
- A new `AlignmentConfig` Pydantic model in
  [src/merxen/config.py](../../src/merxen/config.py).
- `COMPARE` and `VISUALIZE` will receive transformed shape/point layers or
  transforms they can apply on the fly.
- The enriched zarr may be amended with `aligned_*` shape/point layers.

## Why it's deferred

Alignment introduces hyperparameters (landmark choice, elastic vs affine,
tolerance) that interact with segmentation choices. Exposing a locked
comparison pipeline first keeps that decision independent of the rest of
the stack.

## If you need it now

- Compute pair-specific transforms externally.
- Add them to both SpatialData zarrs as `Affine` / `Identity` transforms
  on the shapes and points before running `COMPARE` / `VISUALIZE`.
- Everything downstream will then operate in the shared frame.
