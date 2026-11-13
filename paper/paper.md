---
title: 'SPM 26: Advances in neuroimaging analysis and computational modelling'
tags:
  - neuroimaging
  - neuroscience
  - SPM
  - MRI
  - MEG
  - EEG
  - PET
  - OPM
  - active inference
  - MATLAB
authors:
  - name: Pranay S. Yadav
    orcid: 0000-0003-4692-6506
    affiliation: 1
  - name: Nicholas A. Alexander
    orcid: 0000-0001-5161-684X
    affiliation: 1
  - name: John Ashburner
    orcid: 0000-0001-7605-2518
    affiliation: 1
  - name: Yaël Balbastre
    orcid: 0000-0001-8758-9978
    affiliation: 2
  - name: Gareth Barnes
    orcid: 0000-0002-5396-7712
    affiliation: 1
  - name: Yulia Bezsudnova
    orcid: 0000-0003-2443-7082
    affiliation: 1
  - name: Martina Callaghan
    orcid: 0000-0003-0374-1659
    affiliation: 1
  - name: Barbara Dymerska
    orcid: 0000-0001-6339-2169
    affiliation: 1
  - name: Korbinian Eckstein
    orcid: 0000-0002-4538-7072
    affiliation: 1
  - name: Karl Friston
    orcid: 0000-0001-7984-8909
    affiliation: 1
  - name: Chetan Gohil
    orcid: 0000-0002-0888-1207
    affiliation: 1
  - name: Olivia S. Kowalczyk
    orcid: 0000-0002-6543-290X
    affiliation: "1, 3"
  - name: Vladimir Litvak
    orcid: 0000-0001-8535-7452
    affiliation: 1
  - name: Alberto Mariola
    orcid: 0000-0003-4660-1306
    affiliation: 1
  - name: Johan Medrano
    orcid: 0000-0002-7558-2071
    affiliation: "1, 4"
  - name: George O'Neill
    orcid: 0000-0002-0035-3768
    affiliation: "5, 6"
  - name: Tim M. Tierney
    orcid: 0009-0008-7575-679X
    affiliation: 1
  - name: Peter Zeidman
    orcid: 0000-0003-3610-6619
    affiliation: 1
    corresponding: true
affiliations:
  - name: Department of Imaging Neuroscience, UCL Queen Square Institute of Neurology, University College London, London WC1N 3AR, United Kingdom
    index: 1
    ror: 02jx3x895
  - name: Department of Experimental Psychology, Division of Psychology and Language Sciences, University College London, London WC1H 0AP, United Kingdom
    index: 2
    ror: 02jx3x895
  - name: Department of Neuroimaging, Institute of Psychiatry, Psychology & Neuroscience, King’s College London
    index: 3
    ror: 0220mzb33
  - name: Department of Biological Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA
    index: 4
    ror: 042nb2s44
  - name: Department of Neuroscience, Physiology and Pharmacology, University College London, London, United Kingdom
    index: 5
    ror: 02jx3x895
  - name: FieldLine Medical, Louisville, CO, United States
    index: 6
date: 13 November 2026
bibliography: paper.bib
---

# Summary

Statistical Parametric Mapping (SPM) is an open-source software package for analysing neuroimaging data and testing hypotheses about central nervous system structure and function. Developed over more than 30 years [@zeidman2025spm], `SPM` supports magnetic resonance imaging (MRI), positron emission tomography (PET), electroencephalography (EEG) and magnetoencephalography (MEG), including recordings from optically pumped magnetometers (OPMs). It also provides methods for Bayesian inference and generative modelling of neuronal dynamics and behaviour.

This paper describes `SPM 26`, released in October 2026. Advances since `SPM 25` [@tierney2025spm] include improvements in MRI distortion correction, calibration and interference reduction for OP-MEG, and further developments in M/EEG analysis. Extensions to dynamic causal modelling (DCM) and the dynamic expectation maximisation (DEM) toolbox provide new methods for inference using generative models of brain dynamics, learning and behaviour. We have also expanded numerical regression testing, introduced interoperability with Julia, and improved automated testing of compiled and containerised distributions. Together, these developments extend SPM's analysis and modelling capabilities while helping to maintain the reliability of existing workflows.

# Statement of need

SPM has played a foundational role in the development of statistical methods for cognitive and clinical neuroimaging, including voxel-wise general linear models [@friston1994statistical], topological inference using random field theory [@worsley1996unified], and DCM of interactions between brain regions [@friston2003dynamic]. Building on `SPM 25` [@tierney2025spm], there is a continuing need to extend these capabilities as neuroimaging technologies and computational approaches evolve, while maintaining the reproducibility of existing analyses.

For researchers studying the central nervous system, `SPM 26` addresses three areas of need. First, new neuroimaging technologies and experimental paradigms require improved methods for data preprocessing and analysis. For example, OPMs enable MEG recordings during movement [@boto2018moving; @mellor2023real], creating a need for accurate sensor calibration and interference reduction. Second, developments in computational neuroscience call for more expressive generative models and inference methods to investigate neural dynamics, learning and behaviour. Third, increasingly complex analytical workflows require rigorous numerical validation, reproducible software distributions and interoperability with other programming languages.

# State of the field

Several open-source packages support neuroimaging research, including FSL [@jenkinson2012fsl], AFNI [@cox1996afni] and FreeSurfer [@fischl2012freesurfer] for MRI, and FieldTrip [@oostenveld2011fieldtrip], MNE-Python [@gramfort2013meg], Brainstorm [@tadel2011brainstorm] and EEGLAB [@delorme2004eeglab] for M/EEG. Related packages for computational modelling include TAPAS [@frassle2021tapas] and pymdp [@heins2022pymdp].

SPM differs in its longstanding emphasis on testing hypotheses about the processes underlying observed brain activity through generative modelling and principled statistical inference [@zeidman2025spm]. To this end, it brings together methods for analysing MRI, PET, EEG and MEG data within an extensible framework for Bayesian inference and computational modelling. Continued development within SPM allows new methods to build on these established foundations and become part of existing neuroimaging workflows.

# Software design

SPM is written primarily in MATLAB (81% of source code lines, excluding external dependencies), with computationally intensive operations implemented in C. In `SPM 26`, we have extended its modular system of toolboxes and analysis routines, improved access to scientific software written in other languages, and expanded our testing framework to detect unintended numerical changes between releases.

## Modular extensions

**MRI and statistical inference.** The SCOPE toolbox, introduced in `SPM 25` [@tierney2025spm], corrects geometric distortions in MRI using images acquired with opposite phase-encoding directions [@andersson2003correct]. In `SPM 26`, we extended SCOPE to incorporate an initial displacement map, either supplied by the user or estimated from phase images using ROMEO phase unwrapping [@dymerska2021romeo]. Distortion correction during fMRI preprocessing can also be performed using a static field map, without modelling motion-related changes in the field.

The DICOM conversion utility now supports three-dimensional MR Spectroscopic Imaging (MRSI) data, preserving spatial orientation in five-dimensional NIfTI images. For statistical inference using random-field theory, we added a function to compute Euler characteristic densities for Z, T, F and chi-squared random fields in arbitrary dimensions [@worsley2004unified].

**M/EEG.** `SPM 26` includes a template-based method for removing deep brain stimulation artefacts from EEG evoked responses [@peeters2023electrophysiologic]. For source-space analysis, we added example pipelines to the DAiSS toolbox for studying heartbeat-evoked potentials [@virjee2026heartbeat] using beamformers [@westner2022unified]. Bayesian spectral decomposition (BSD), previously described in the `SPM 25` paper, is now included as a toolbox [@medrano2025bsd]. BSD estimates parametric models of neural power spectra and uses Bayesian model comparison to assess whether particular oscillatory peaks are supported by the data.

**OP-MEG.** Accurate analysis of OPM recordings depends on knowledge of sensor positions, orientations and gains, as well as effective suppression of environmental interference. In `SPM 26`, we introduced coil-based calibration, in which known magnetic fields generated by external coils are used to estimate sensor geometry and sensitivity [@hill2025determining]. We also added an implementation of signal-space projection [@uusitalo1997signal] that supports arrays with unequal noise variance across channels, as can occur in multi-axis OPM systems. Anatomically derived two-dimensional sensor layouts have also been extended to represent tangential channels [@alexander2025anatomically].

**Dynamic causal modelling.** Conventional DCMs describe effective connectivity through parameters that are typically assumed to remain constant over the period being modelled. `SPM 26` adds demonstrations of time-varying effective connectivity for event-related potentials and cross-spectral densities, using temporal basis functions to represent slow changes in synaptic efficacy [@medrano2024dynamic]. For DCM of event-related potentials, users can now select from SPM's built-in exponential integrators or provide their own integration routines.

**Active inference and behavioural modelling.** We substantially expanded the DEM toolbox [@friston2008dem] with new methods and demonstrations for generative modelling, learning, inference and decision-making. New renormalising generative models and structure-learning methods explore how hierarchical representations can be learned from sensory data, with examples involving images, video and sound [@friston2025pixels]. New methods allow agents to consider sequences of future actions and their expected outcomes when selecting policies for planning problems such as maze navigation and the Tower of Hanoi [@friston2021sophisticated]. Other developments include non-equilibrium steady-state models and methods based on Helmholtz–Hodge decompositions for studying dynamical systems, forecasting and self-organisation [@friston2026free]. We also revised core active inference routines [@friston2017active], including policy representations, hierarchical model specification and Bayesian model reduction [@friston2017curiosity].

## Interoperability and reuse

We have extended SPM's interfaces to other programming languages, allowing researchers to use specialised and computationally optimised methods developed outside MATLAB. In `SPM 26`, we introduced `spm_julia`, which enables MATLAB-based workflows to execute code written in Julia [@bezanson2017julia]. The interface manages a local Julia installation and its package dependencies, providing access to Julia methods without requiring users to configure a separate environment. Its first application is ROMEO [@dymerska2021romeo], a fast and accurate phase-unwrapping method used for field-map estimation in SCOPE, suitable for large-scale fMRI studies. The `spm-python` interface, introduced in `SPM 25` [@tierney2025spm], is now available through PyPI, allowing SPM to be used from Python-based workflows.

SPM has long used components of FieldTrip for M/EEG data input and forward modelling, drawing on established implementations maintained by the wider research community. In `SPM 26`, we introduced a continuous integration workflow that tracks FieldTrip's upstream release branch and incorporates updates once tests succeed. This allows SPM to benefit from ongoing improvements to FieldTrip while checking that updates remain compatible with existing SPM functionality.

## Numerical validation and reproducibility

Automated testing has expanded in `SPM 26`, with the number of test cases increasing from 211 to 337. Because changes to numerical routines can alter model predictions and parameter estimates without causing software failures, new regression tests compare numerical outputs against stored reference values to detect unintended changes.

A new DCM test suite includes parameterised tests spanning nine neural mass models and three numerical integration schemes. These tests check each model's stability and steady-state initialisation. They also assess whether models reproduce simulated responses and whether inversion of models from fixed datasets yields estimates that agree with stored reference values within specified tolerances.

Additional regression tests evaluate M/EEG analysis workflows, including a pipeline using the Wakeman and Henson dataset [@wakeman2015multi] that covers preprocessing through source reconstruction. Automated workflows run tests across supported configurations of Linux, macOS and Windows, using MATLAB releases from R2020b onwards. Changes to numerical behaviour and processing defaults that may affect reproducibility are documented in the release notes.

## Build and distribution

As in `SPM 25`, we test and distribute `SPM 26` in three forms, each supporting GUI and command-line use. The standard MATLAB package runs on personal computers and high-performance computing (HPC) systems but requires a MATLAB licence. Standalone packages provide compiled versions of SPM for Linux, Windows and macOS, including Intel and Apple Silicon architectures, without requiring a MATLAB licence. Container images for Docker, Podman and Singularity provide a further licence-free option for deployment across computing environments. Maintaining these distributions increases build and testing requirements but makes SPM accessible to users with different software and computing resources.

Development takes place publicly on GitHub. Between releases `25.01.02` and `26.10.01`, 24 people contributed about 350 commits, including seven from outside the SPM organisation. SPM is distributed under the GNU General Public License (version 2 or later), with documentation available at https://www.fil.ion.ucl.ac.uk/spm/docs/.

# Research impact statement

SPM is widely used in cognitive, clinical and computational neuroscience. To assess the research community served by the software, we examined publications reporting use of `SPM 25` since its release in January 2025. We identified 76 publications and software records, including 41 journal articles, 67 of them with no authors in common with the SPM development team or the `SPM 25` paper (full list available in the [research impact bibliography](https://github.com/spm/spm/blob/spm-release-paper/paper/research_impact_bibliography.md)). These studies illustrate the breadth of SPM's applications, spanning task fMRI [@knobloch2026resolving; @memis2026object], resting-state fMRI [@chuda2026early], structural MRI [@jiskoot2026automated], PET [@chiaravalloti2026plasma; @li2026total], DCM [@taniguchi2026task; @kvistad2026effects] and functional near-infrared spectroscopy [@levy2026shedding].

`SPM 26` builds on this established research use while extending support for emerging applications. In naturalistic neuroimaging, its OP-MEG developments support studies of brain activity during movement and in participants for whom conventional MEG may be difficult [@hill2019tool; @vivekananda2020optically]. At the same time, its developments in DCM and the DEM toolbox aim to improve inference using generative models across a range of applications.

# AI usage disclosure

We used generative AI in two areas. First, Claude (Anthropic; models Opus 4.7, 4.8, 5 and 5.5) assisted with writing code and commit messages for the build, testing and continuous integration infrastructure. All AI-assisted code was reviewed and verified by an SPM developer. Second, we used Claude (Opus 5.5) to help draft the `SPM 26` release notes, which were extensively reviewed, edited and verified by the SPM development team before publication.

# Conflict of interest

GCO is employed by FieldLine Medical (FLM), a company commercialising OP-MEG technology. FLM had no editorial control over the content of this publication.

# Acknowledgements

A comprehensive list of `SPM` authors is in the file `AUTHORS.txt` supplied with the software. We are also grateful to the IT Team at the UCL Department of Imaging Neuroscience for their ongoing support.

Yael Balbastre is funded by a Royal Society Newton International Fellowship (NIF\\R1\\232460). Olivia S. Kowalczyk is supported by the King's Prize Fellowship. Nicholas A. Alexander, Barbara Dymerska, Alberto Mariola and Johan Medrano are supported by the Discovery Research Platform for Naturalistic Neuroimaging funded by Wellcome \[226793/Z/22/Z\]. Barbara Dymerska is also partially funded by the Medical Research Council \[UKRI226\]. Tim M. Tierney is funded by an Epilepsy Research UK fellowship (FY2101). Peter Zeidman is funded by an MRC Career Development Award \[MR/X020274/1\].

# References
