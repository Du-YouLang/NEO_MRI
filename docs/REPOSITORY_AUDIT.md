# Repository audit

## Scope

This review starts from remote `main` at `1ec846a` and the five published feature-branch tips listed in the [branch index](../README.md#branches-and-validation). It follows each public entry point through configuration, geometry, reporting and tests. The existing [Ubuntu CI run](https://github.com/Du-YouLang/NEO_MRI/actions/runs/43363179281) applies to the baseline commit. Windows tests and the synthetic tutorials exercise this review's changes. The review did not run a new patient scan, native MRtrix3 tracking, a real HCP import or a production fMRIPrep container.

## Findings

| Priority | Location | Finding and status |
| --- | --- | --- |
| High | [BOLD saved-map operations](../bold/bold_stage.py) | Two equal-sized arrays can represent different physical locations. The code now checks that each saved statistical map is 3D and on the mask's physical grid before comparison or fixed-effects calculation. A shifted-affine regression exercises both paths. |
| Medium | [T1 planning](../code/t1_preprocess.py) | The planner could hash an out-of-BIDS T1 before the launcher rejected it. Planning now validates the selected file's path and dimensionality before hashing it. A regression exercises an out-of-tree copy with the same content. |
| Medium | [Jupyter subprocesses](../notebooks/Task_BOLD.ipynb) | A notebook may export a Matplotlib inline backend unavailable to the subprocess that renders HTML. The staged BOLD and tractography notebooks now set a headless subprocess backend and UTF-8 explicitly. The BOLD event override cell also exposes the required reviewer and reason fields. |
| High | [Tractography atlas mapping](../tractography/geometry.py) | The MNI6 bridge and archived anatomy H5 are checked by hash before use. Their direction and anatomical correspondence still require independent landmarks or approved point-transform evidence. The CST work found a template-bridge pair with a roughly 5.7 mm round-trip 95th-percentile error in its earlier local validation. That result applies to the tested pair and points. This audit did not change the bridge or bless the direct pull. |
| High | [Lesion load archived input chain](../lesionload/spatial.py) | `archived_explicit_chain` records the review file and its hash. It does not parse the verdict or bind it semantically to the visit and the supplied transforms. A reviewer must verify those links. The manifest-based structural path provides more evidence. |
| Medium | [Native MRtrix execution](../tractography/workflow.py) | Completion parsing depends on MRtrix3 progress text and TCK headers. The README distinguishes synthetic and saved-streamline checks from qualification of the native launcher in a real environment. |
| Medium | [BOLD statistics](BOLD_AB.md) | The dynamic search and bootstrap are implemented and covered by generated examples. A and B have different inferential scope. B is exploratory, and A's 1,214 multiplier is the project's convention. Clinical or population calibration requires separate evidence. |
| Maintenance | [CI configuration](../.github/workflows/tests.yml) | The reported Node.js action warning and Ubuntu runner notice are separate from pytest failures. The latest baseline CPU workflow succeeded. Runner and action-version changes are outside this audit. |

The BOLD fit fingerprint includes the full `bold_stage.py` source. This change creates a new fit identity under that policy. Use the original commit to inspect earlier cached fits. Existing outputs are preserved.

## Local acceptance

- Windows Python 3.12, 56 pytest tests passed. The structural wheel was built and its 22 packaged modules were verified.
- All five notebooks executed with synthetic data or the default inactive production controls. The three tutorials produced HTML. The optional OpenNeuro exercise verified the exact TSV hash and read 185 event rows.
- The 12 bilingual README flowcharts were parsed and rendered to SVG in a DOM test. The generated native-mask, registration blink, multimodal and lesion-load HTML controls passed their Node DOM/canvas simulations. Those tests exercise logic, while a reviewer evaluates actual alignment and display in the browser.

## Public-data and rights review

The six remote branches and their 16 reachable commits were checked for the private case identifiers and machine-local paths, Git-tracked neuroimaging, archives and stored notebook outputs. No private scans or matching identifiers were found. A small set of internal case-derived summary counts was present in the BOLD guide. The current documentation uses synthetic exercises and source-backed public examples. Published historical commits retain their original contents.

The CST and DoubleScanner maintainers independently checked the sources, public-versus-local method boundaries and validation status. [Acknowledgments](ACKNOWLEDGMENTS.md) distinguishes dependencies, external assets, and related local experiments. The team's publication licence remains pending as recorded in [NOTICE](../NOTICE.md).
