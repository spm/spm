---
title: 'SPM 25: open source neuroimaging analysis software'
tags:
  - neuroimaging
  - neuroscience
  - SPM
  - MRI
  - MEG
  - EEG
  - PET
authors:
  - name: John Ashburner
    affiliation: 1	
  - name: Nicholas Alexander
    affiliation: 1	
  - name: Nicole Labra Avila
    affiliation: 1	
  - name: Yaël Balbastre
    affiliation: 1	
  - name: Gareth Barnes
    affiliation: 1	
  - name: Yulia Bezsudnova
    affiliation: 1	
  - name: Mikael Brudfors
    affiliation: 1	
  - name: Korbinian Eckstein
    affiliation: 2	
  - name: Guillaume Flandin
    affiliation: 1	
  - name: Karl Friston
    affiliation: 1	
  - name: Amirhossein Jafarian
    affiliation: 1		
  - name: Olivia Kowalczyk
    affiliation: 1	
  - name: Vladimir Litvak
    affiliation: 1	
  - name: Johan Medrano
    affiliation: 1	
  - name: Stephanie Mellor
    affiliation: 1	
  - name: George O'Neill
    affiliation: 1	
  - name: Thomas Parr
    affiliation: 1	
  - name: Adeel Razi
    affiliation: 1	
  - name: Tim Tierney
    affiliation: 1	
  - name: Ryan Timms
    affiliation: 1	
  - name: Peter Zeidman
    orcid: 0000-0003-3610-6619
    corresponding: true
    affiliation: 1
affiliations:
 - name: University College London, UK
   index: 1
   ror: 02jx3x895
 - name: Independent Researcher, Germany
   index: 2
date: 18 December 2024
bibliography: paper.bib

---

# Summary

Statistical Parametric Mapping (SPM) is an integrated set of methods for testing hypotheses about the brain's structure and function, using data from medical imaging devices. The open source software package implementing these methods, `SPM`, has been in continuous development for more than 30 years by an international community of developers, working on the core software and on toolboxes that extend its functionality. This paper reports the release of `SPM 25`, a major new version of the software that incorporates novel analysis methods, optimisations of existing methods, as well as improved practices for open science and software development.

# Statement of need

`SPM` introduced many of the statistical methods that underpin cognitive and clinical neuroimaging research today, including:

- The voxel-wise application of General Linear Models (GLMs) to imaging data
- Correction for multiple comparisons across space using topological methods (Random Field Theory, RFT).
- Convolution modelling of functional MRI (fMRI) signals
- Event-related fMRI for cognitive neuroscience
- Voxel-Based morphometry (VBM) for detecting changes in anatomy
- Dynamic Causal Modelling (DCM) for state-space modelling using variational Bayesian methods
- Source localisation for M/EEG data using variational Bayesian methods

What unifies these methods are a series of core principles: the use of generative models that explain (rather than simply describe) data, the well-motivated application of parametric statistics for testing hypotheses, and open science practices that have meant that `SPM` has been freely available since its inception.

This major new release of `SPM` addresses a series of needs in the neuroimaging community, set out below.

## Open development

`SPM` was previously developed and tested using a private Subversion server within University College London. To enable community engagement in the future of `SPM` and to increase transparency, development has now moved to a public Github repository. `SPM 25` is the first release following this change. The key advantages so far have been:

- Introducing automated unit and regression tests across platforms.
- Issue tracking
- Automating the build process to easily generate and release source code and compiled versions.

## Documentation and training

`SPM 25` is accompanied by a new [documentation website]((https://www.fil.ion.ucl.ac.uk/spm/docs/)) with step-by-step tutorials on all of SPM's main features, as well as freely available video recordings of lectures from previous SPM courses covering the mathematical theory. This website is under development and is intended to unify the disparate sources of documentation that were previously available, such as the PDF manual.

## Major new features 

`SPM 25` includes 10 years of new developments since the last major release (`SPM 12`, dated 1st October 2014). This section highlights some of the most significant new features.

### MRI

New features for fMRI and multimodal data include:

- Parametric Empirical Bayes (PEB) for DCM [@friston2016bayesian]. Introduces random effects modelling of neural connectivity parameters, enabling people to test hypotheses about the similarities and differences among research participants.

- Multi-Brain Toolbox [@brudfors2020flexible]. Generates population average-shaped brains, enabling more precise spatial normalisation with the option to automatically label brain structures.

- SCOPE Toolbox. Generates MRI fieldmaps using phase-encode-reversed pairs of MRI images (blip-up and blip-down images), similar to the Topup tool in FSL.

### M/EEG

Some new features specific to M/EEG data are:

- Methods for spectral decomposition - `SPM 25` offers both an implementation of an existing approach called FOOOF (Specparam) in the MEEGtools toolbox, based on code from Brainstorm [@donoghue2020parameterizing], as well as a new Bayesian implementation that introduces formal statistical testing, called Bayesian Spectral Decomposition (BSD) [@medrano2024bsd].

- Support for fusion of different MEG sensor type and EEG in beamforming with pre-whitening [@westner2022unified].

- Support for MEG BIDS for specification of events, channels and fiducials [@westner2022unified].

### OPMs

MEG using Optically Pumped Magnetometers (OPMs) enables free movement of the head and body during MEG recordings [@boto2018moving]. This opens MEG to new scientific applications (e.g., experiments involving free movement) and new study populations (e.g., children). Developing the analysis tools for OPM data is a major focus for SPM, with recently added features including:

- File IO for all major OPM manufacturers (Quspin, Cerca, Mag4Health, Fieldline).
- Methods to simulate arbitrary OPM arrays of differing densities and vector measurements.
- OPM interference cancellation algorithms for low channel systems: Homogeneous Field Correction [@tierney2021modelling].
- OPM interference cancellation algorithms for large channel systems: Adaptive Multipole Models [@tierney2024adaptive].

### Making SPM accessible without MATLAB

Approximately 90% of `SPM 25` is written in MATLAB and the remainder is C++. This code has been highly optimised and thoroughly tested over its 30 years of development. We have carefully considered how to capitalise on the stability of the `SPM` software, while introducing accessibility for people who do not have access to a MATLAB license, or who prefer to write their analysis code in other languages. 

Our strategy is as follows:

- `SPM 25` will be the first version of SPM to be fully accessible from the Python programming language, without requiring MATLAB, using a new Python wrapper called [spm-python](https://github.com/spm/spm-python). This is in the final stages of testing and will be released in the first quarter of 2025.

- [SPM Standalone](https://www.fil.ion.ucl.ac.uk/spm/docs/installation/standalone/) is the compiled version of SPM that can be run from the command line without a MATLAB license. It is now generated automatically with each new release, as part of SPM's Github-based build process. This enables people to run neuroimaging analyses from command line scripts.

- [Docker and Singularity containers](https://www.fil.ion.ucl.ac.uk/spm/docs/installation/containers/) are additionally provided and are now generated automatically as part of SPM's Github build process.

# Boilerplate text to delete

## Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

## Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

## Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

JM is supported by the Discovery Research Platform for Naturalistic Neuroimaging funded by Wellcome
[226793/Z/22/Z]. PZ is funded by an MRC Career Development Award [MR/X020274/1].

# References