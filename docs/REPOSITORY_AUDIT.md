# Repository audit

## Scope

This review starts from `main` at `1ec846a` and the five feature-branch tips listed in the [branch index](../README.md#branches-and-validation). It covers public entry points, configuration, geometry, reporting and tests. The [baseline Ubuntu run](https://github.com/Du-YouLang/NEO_MRI/actions/runs/43363179281) and [initial spatial-review PR checks](https://github.com/Du-YouLang/NEO_MRI/pull/6/checks) apply to their tested revisions. GitHub Actions runs again on the integrated revision after it is pushed.

## Findings

| Priority | Location | Finding and status |
| --- | --- | --- |
| High | [BOLD saved-map operations](../bold/bold_stage.py) | Two equal-sized arrays can represent different physical locations. The code now checks that each saved statistical map is 3D and on the mask's physical grid before comparison or fixed-effects calculation. A shifted-affine regression exercises both paths. |
| Medium | [T1 planning](../code/t1_preprocess.py) | The planner could hash an out-of-BIDS T1 before the launcher rejected it. Planning now validates the selected file's path and dimensionality before hashing it. A regression exercises an out-of-tree copy with the same content. |
| Medium | [Jupyter subprocesses](../notebooks/Task_BOLD.ipynb) | A notebook may export a Matplotlib inline backend unavailable to the subprocess that renders HTML. The staged BOLD and tractography notebooks now set a headless subprocess backend and UTF-8 explicitly. The BOLD event override cell also exposes the required reviewer and reason fields. |
| High | [Tractography atlas mapping](../tractography/geometry.py) | Gen6 mapping retains the pinned point chain via MNI2009c. Reviewed mapping and real streamline screening now bind transform directions, reference grids, case landmarks and a scope-specific approval to the inputs. A prior round-trip p95 of about 5.7 mm applies to the tested pair of independently published transforms and the sampled SMATT domain. Each new visit requires anatomical review. |
| High | [Lesion load inputs](../lesionload/spatial.py) | New runs verify an approval that binds the actual visit, original lesion, images, transforms and, where applicable, the frozen T2 branch. Optional FA requires its own spatial and model evidence. Existing approvals can be reused within their documented scope. |
| Medium | [Native MRtrix execution](../tractography/workflow.py) | Completion parsing depends on MRtrix3 progress text and TCK headers. The README distinguishes synthetic and saved-streamline checks from qualification of the native launcher in a real environment. |
| Medium | [BOLD statistics](BOLD_AB.md) | The dynamic search and bootstrap are implemented and covered by generated examples. A and B have different inferential scope. B is exploratory, and A's 1,214 multiplier is the project's convention. Clinical or population calibration requires separate evidence. |
| Maintenance | [CI configuration](../.github/workflows/tests.yml) | The reported Node.js action warning and Ubuntu runner notice are separate from pytest failures. The latest baseline CPU workflow succeeded. Runner and action-version changes are outside this audit. |

The BOLD fit fingerprint includes the full `bold_stage.py` source. This change creates a new fit identity under that policy. Use the original commit to inspect earlier cached fits. Existing outputs are preserved.

## Local acceptance

- The integrated checkout passed all 63 tests in the Windows Python 3.12 CPU suite. The structural wheel contained the expected 22 modules, and the installed structural CLI and T2 policy worked from outside the checkout. All five notebooks executed with synthetic data or the default inactive production controls. The three tutorials produced HTML. The optional OpenNeuro exercise verified the exact TSV hash and read 185 event rows.
- All 12 bilingual README flowcharts were parsed and rendered to SVG in a DOM test. Node DOM/canvas checks exercised the generated native-mask, registration blink, multimodal and lesion-load viewers. The lesion-load test now checks values by tract name across the three methods, whose alphabetical display orders differ.
- The spatial revision reproduced 12 existing brainstem labels across two archived mappings. Three archived lesion-load cases reproduced 18 lesion and field-of-view maps and 22,062 table rows. The replay used the existing images and transforms and provides reproducibility evidence. A native MRtrix3 run, HCP imaging import and production fMRIPrep container run remain outside this acceptance.

## Public-data and rights review

The baseline repository, published feature branches and pending release were checked for private-case identifiers, machine-local paths, Git-tracked neuroimaging, archives and stored notebook outputs. The scan found no private scans or matching identifiers. The current documentation uses synthetic exercises and source-backed public examples. Published historical commits retain their original contents.

The CST and DoubleScanner maintainers independently checked the sources, public-versus-local method boundaries and validation status. [Acknowledgments](ACKNOWLEDGMENTS.md) distinguishes dependencies, external assets, and related local experiments. The team's publication licence remains pending as recorded in [NOTICE](../NOTICE.md).
