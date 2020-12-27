# SpatialAttentionNetwork
A network model of selective attention in the auditory cortex, implemented in the DynaSim framework.

This projects simulates 3 types of phenomena:
* Selective attention in the spatial domain (Lee & Middlebrooks, 2011)
* Selective attention in the spectral domain (Fritz et al., 2003, and Atiani et al., 2009)
* Functional implications of selective attention in the spatial and spectral domains

## Spatial Domain Simulations
The script `runSpatialAtten_main.m` runs simulations offering an explanation to the observed attention-induced changes in receptive fields in the cat A1 (Lee & Middlebrooks 2003).
We offer three types of mechanisms to explain these observations
* changes due to the cholinergic system
* changes due to top-down attention
* combination of the two

Comment and uncomment the appropriate line to load the parameters for each simulation.
See description in the script for more information.

## Spectral Domain Simulations
The script `runFreqAttenNetwork.m` runs simulations to offer explanations to sme of the observed attention-induced changes in receptive fields.

We simulated the following:
* Sharpening of receptive fields (Fig 3a, Atiani et al., 2009)
* Weakening of receptive fields (Fig 3b, Atiani et al., 2009)
* Emergence of new receptive area (Fritz et al., 2003)

Comment and uncomment the appropriate line to load the parameters for each simulation.
See description in the script for more information.

## Functional Simulations
We also simulated functional uses for the AIM network:
* `runSpatialAttenApplication_main.m` simulates spatial attention
* `runHarmonicAttenNetwork_main.m` simulates attending to harmonics of a speaker's f0
