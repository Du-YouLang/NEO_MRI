# Public data and tutorials

The notebooks ship without human images. Start with the generated BOLD, tractography and lesion load examples, then use one of the sources below to practice importing real data. Keep downloads, converted data, configuration copies and reports outside the Git checkout.

| What to learn | Dataset | Files to choose |
| --- | --- | --- |
| T1 and T2 review | [HCP Young Adult, 2025 release](https://humanconnectome.org/study/hcp-young-adult/document/hcp-young-adult-2025-release) | 3T T1w and T2w for the same participant. Choose unprocessed inputs to test your own extraction and registration. |
| DWI, CST and CRP | [HCP Young Adult, 2025 release](https://humanconnectome.org/study/hcp-young-adult/document/hcp-young-adult-2025-release) | A corrected 3T diffusion volume and its matching bvals and bvecs from the same series, plus the native T1 and a reviewed DWI–T1 transform when mapping atlas gates. |
| Motor task fMRI | [HCP Young Adult task protocol](https://www.humanconnectome.org/hcp-protocols-ya-task-fmri) | The MOTOR LR and RL acquisitions with their own EV files and metadata. These are movement-execution runs. |
| BIDS task fMRI and timing | [OpenNeuro ds000117 v2.0.0](https://doi.org/10.18112/openneuro.ds000117.v2.0.0) | `sub-01/ses-mri/task-facerecognition/run-01` BOLD, JSON and events, plus same-session T1 and applicable distortion-correction inputs. |
| A genuine stroke lesion | [ARC / OpenNeuro ds004884 v1.0.2](https://doi.org/10.18112/openneuro.ds004884.v1.0.2) | A same-session T2 and expert lesion mask, such as `sub-M2026`. Check the T2 acquisition and voxel geometry against the actual annotation. |
| 7 T diffusion input practice | [MICA-PNI derivatives](https://osf.io/download/ksnht/?version=1) | The release's documented single-participant diffusion derivatives, with provenance retained. See the [dataset paper](https://doi.org/10.1038/s41597-025-04863-7). |

## HCP Young Adult

Use [ConnectomeDB and the HCP data-use terms](https://www.humanconnectome.org/study/hcp-young-adult/data-use-terms) to obtain the appropriate packages. Access and redistribution conditions apply. Record the release, package, participant, input filenames and acquisition identity. The [release package contents](https://humanconnectome.org/storage/app/media/documentation/HCP-YA2025/HCP-YA_2025_Release_Appendix_III.pdf) document the native and standard-space files.

1. **T1/T2.** Identify the original T1w and T2w for one participant. Record the phase-encoding-independent structural geometry, voxel size and scan identifiers. Use the [T1 guide](T1_MASK_REVIEW.md) to prepare the full-head N4 T1 and review the brain mask, then point the [T2 configuration](../configs/examples/t2_v1.example.json) at the original T2 and the approved T1 chain.
2. **Diffusion.** In the 3T recommended package, use `T1w/Diffusion/data.nii.gz`, `bvals`, `bvecs` and `nodif_brain_mask.nii.gz`. Check the gradient count and order against the volume count. The HCP output is already aligned to the participant's T1w space. Verify which T1 it matches before supplying a transform and retain the HCP gradient nonlinearity information if the fitting protocol requires it. The [diffusion tutorial](../notebooks/Tractography_Tutorial.ipynb) shows how to provide the DWI, masks and transform. For the current native single-tissue route, select a reviewed b-value shell and its matching gradient rows as described in the [tractography guide](TRACTOGRAPHY.md).
3. **Task fMRI.** Open the EV files for one MOTOR acquisition at a time. The three numeric columns are onset in seconds, duration in seconds, and weight. Preserve the weight in your intake manifest and review how to represent it in your chosen model. Stage a versioned `events.tsv` with `onset`, `duration` and `trial_type`, using the EV file identity for the condition. Keep LR and RL separate and verify each run's reference clock, image and metadata. Check that the event and contrast specification matches the scientific estimand before considering a joint fit.
4. **BOLD derivatives.** The current [BOLD intake](BOLD_AB.md) expects registered fMRIPrep derivatives, confounds and preprocessing evidence. HCP's supplied processing follows a different contract. Start from a reviewed BIDS conversion and run fMRIPrep externally to use the existing intake. Document the conversion, retain the source HCP files, and verify phase-encoding metadata and run timing against the release protocol.

HCP is a healthy reference cohort. Use an explicitly synthetic mask or the ARC example below to exercise lesion handling. This repository does not distribute HCP images or provide an HCP conversion pipeline.

## OpenNeuro BIDS task exercise

Start with the [pinned dataset README](https://github.com/OpenNeuroDatasets/ds000117/blob/2.0.0/README). The example run has 208 saved volumes at TR 2 s. The publisher removed two dummy volumes before the saved time origin. Preserve the original events and check both the raw clock and the actual derivative timing metadata. With slice-timing correction, the code uses the derivative's `StartTime` together with TR. Avoid adding the dummy-volume offset a second time.

Use a reviewed, documented event selection for static GLM and FIR. The released events include stimulus and response rows. For DoubleScanner, define and review the specific dynamic condition and scan timing first. The BOLD tutorial exercises that interface using generated tasks. The public example is available for real-data intake practice and was not rerun by this repository audit.

## ARC lesion exercise

Use the dataset documentation and the `derivatives/lesion_masks` tree to locate the matching T2 and binary lesion annotation. Confirm affine, shape and visit identity, then review the original label against T2 slices. The [lesion interface](LESION.md) exports on its configured reference grid. Review a mapped copy if using it on a T1 grid. For [lesion load](LESIONLOAD.md), keep the native T2 mask and provide the actual reviewed transforms and atlas versions. In ARC, a task event table can cover a different recording duration from a sparse BOLD acquisition. Check both the scan length and the event order before fitting any task.

The [acknowledgments](ACKNOWLEDGMENTS.md) give the source citations. Use the HCP provider's required acknowledgment when reporting results from its data.
