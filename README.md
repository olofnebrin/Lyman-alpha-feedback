## Repository for project on Lyman-α radiation pressure feedback
<div style="text-align: center;">
    <img src="https://github.com/user-attachments/assets/fa23c18b-6641-4bc1-beed-804950a428bb" alt="Schematic Overview Lya" width="550"/>
</div>

### Purpose of repository

This repository contains files and data related to the following paper on Lyman-α (Lyα) feedback (Nebrin et al. 2024 for short below):

[Nebrin et al. (2024), 'Strong Lyman-α feedback at Cosmic Dawn: Implications for the first stars and galaxies, star clusters, and direct-collapse black holes'](INSERT).

Authors & contributors: Olof Nebrin, Aaron Smith, Kevin Lorinc, Johan Hörnquist, and Åsa Larsson. 

### Description of folders

#### M_F_fit:

This folder contains a Python implementation of the fit to the Lyα force multiplier, as given the paper, as well as related quantities 
(e.g. the destruction probability). An example plot is included too. If you use these fits or results in your work, please cite [Nebrin et al. (2024)](INSERT).

#### MCRT_forcemultiplier:

This folder contains the Monte Carlo radiative transfer (MCRT) data for the Lyα force multiplier, as plotted in Fig. 5 in Nebrin+ (2024). Each data file is named in a self-explanatory manner with respect to Fig. 5 in Nebrin+ (2024). In each file, the first column is a_v*tau_cl, and the second column the force multiplier. If you use this data, cite [Nebrin et al. (2024)](INSERT), and [Smith et al. (2015)](https://ui.adsabs.harvard.edu/abs/2015MNRAS.449.4336S/abstract). 

#### H_H_cross_section:

The file in this folder contains the data for the cross-section for H(2s) + H(1s) -> H(2p) + H(1s). The first column is collision energy (in eV), and the second column the cross-section (in cm^2). If you use this data, cite [Nebrin et al. (2024)](INSERT), [Hörnquist et al. (2022)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.106.062821), and [Hörnquist et al. (2023)](https://journals.aps.org/pra/abstract/10.1103/PhysRevA.108.052811).
