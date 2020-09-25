# Multi-Brain: Unified segmentation of population neuroimaging data

## Overview

This repository contains the Multi-Brain (MB) model, which has the general aim of integrating a number of disparate image analysis components within a single unified generative modelling framework (segmentation, nonlinear registration, image translation, etc.). The model is described in Brudfors et al [2020], and it builds on a number of previous works. Its objective is to achieve diffeomorphic alignment of a wide variaty of medical image modalities into a common anatomical space. This involves the ability to construct a "tissue probability template" from a population of scans through group-wise alignment [Ashburner & Friston, 2009; Blaiotta et al, 2018]. Diffeomorphic deformations are computed within a *geodesic shooting* framework [Ashburner & Friston, 2011], which is optimised with a Gauss-Newton strategy that uses a multi-grid approach to solve the system of linear equations [Ashburner, 2007]. Variability among image contrasts is modelled using a much more sophisticated version of the Gaussian mixture model with bias correction framework originally proposed by Ashburner & Friston [2005], and which has been extended to account for known variability of the intensity distributions of different tissues [Blaiotta et al, 2018]. This model has been shown to provide a good model of the intensity distributions of different imaging modalities [Brudfors et al, 2019]. Time permitting, additional registration accuracy through the use of shape variability priors [Balbastre et al, 2018] will also be incorporated.

## Dependencies
The algorithm is developed using MATLAB and relies on external functionality from the SPM12 software:
* **SPM12:** Download from https://www.fil.ion.ucl.ac.uk/spm/software/spm12/.
* **Shoot toolbox:** Add Shoot folder from the toolbox directory of the SPM source code.
* **Longitudinal toolbox:** Add Longitudinal folder from the toolbox directory of the SPM source code.

## References
* Ashburner J, Friston KJ. Unified segmentation. Neuroimage. 2005 Jul 1;26(3):839-51.
* Ashburner J. A fast diffeomorphic image registration algorithm. Neuroimage. 2007 Oct 15;38(1):95-113.
* Ashburner J, Friston KJ. Computing average shaped tissue probability templates. Neuroimage. 2009 Apr 1;45(2):333-41.
* Ashburner J, Friston KJ. Diffeomorphic registration using geodesic shooting and Gauss–Newton optimisation. NeuroImage. 2011 Apr 1;55(3):954-67.
* Blaiotta C, Cardoso MJ, Ashburner J. Variational inference for medical image segmentation. Computer Vision and Image Understanding. 2016 Oct 1;151:14-28.
* Blaiotta C, Freund P, Cardoso MJ, Ashburner J. Generative diffeomorphic modelling of large MRI data sets for probabilistic template construction. NeuroImage. 2018 Feb 1;166:117-34.
* Balbastre Y, Brudfors M, Bronik K, Ashburner J. Diffeomorphic brain shape modelling using Gauss-Newton optimisation. In International Conference on Medical Image Computing and Computer-Assisted Intervention 2018 Sep 16 (pp. 862-870). Springer, Cham.
* Brudfors M, Ashburner J, Nachev P, Balbastre Y. Empirical Bayesian Mixture Models for Medical Image Translation. In International Workshop on Simulation and Synthesis in Medical Imaging 2019 Oct 13 (pp. 1-12). Springer, Cham.
* Brudfors M, Balbastre Y, Flandin G, Nachev P, Ashburner J. Flexible Bayesian Modelling for Nonlinear Image Registration. To appear in: International Conference on Medical Image Computing and Computer-Assisted Intervention 2020

## Acknowledgements
This work was funded by the EU Human Brain Project’s Grant Agreement No 785907 (SGA2).
