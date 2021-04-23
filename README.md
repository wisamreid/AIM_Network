# AIM Network - Attentional inhibition modulation network.
A network model of selective attention in the auditory cortex, implemented using the [DynaSim framework](https://github.com/DynaSim/DynaSim).
This repo accompanies our [preprint on bioRxiv](https://www.biorxiv.org/content/10.1101/2020.12.10.419762v1.abstract).

This project simulates the tuning curves changes of the A1 due to selective attention in several scenarios:
* Selective attention in the spatial domain (Lee & Middlebrooks, 2011)
* Selective attention in the spectral domain (Fritz et al., 2003, and Atiani et al., 2009)
* Functional implications of selective attention in the spatial and spectral domains

We used either pure-tones or white gaussian noise to drive the network. Inputs to the simulations are simulated spike trains generated using the [BOSSA algorithm](https://github.com/kfchou/BOSSA).

Functional applications of the AIM network uses speech stimuli from the [CRM corpus](https://github.com/LABSN/expyfun-data/tree/master/crm). In the spatial attention examples, speech stimuli were spatialized with [KEMAR head-related transfer functions](https://github.com/kfchou/BOSSA/tree/master/HRTF_40k).

## Spatial Domain Simulations
The script `runSpatialAtten_main.m` runs simulations offering an explanation to the observed attention-induced changes in receptive fields in the cat A1 (Lee & Middlebrooks 2003).
We offer three types of mechanisms to explain these observations
* changes due to the cholinergic system
* changes due to top-down attention
* combination of the two

Comment and uncomment the appropriate line to load the parameters for each simulation.
See description in the script for more information.

## Spectral Domain Simulations
The script `runFreqAttenNetwork_main.m` runs simulations to offer explanations to sme of the observed attention-induced changes in receptive fields.

We simulated the following:
* Sharpening of receptive fields (Fig 3a, Atiani et al., 2009)
* Weakening of receptive fields (Fig 3b, Atiani et al., 2009)
* Emergence of new receptive area (Fig 2a, Fritz et al., 2003)

Comment and uncomment the appropriate line to load the parameters for each simulation.
See description in the script for more information.

## Functional Simulations
We also simulated functional uses for the AIM network:
* `runSpatialAttenApplication_main.m` simulates attending to a specific location in space
* `runHarmonicAttenNetwork_main.m` simulates attending to harmonics of a speaker's f0
