# EicC Sullivan Process Reproduction

## Objective
This repository contains the simulation framework, physics validation, and plotting scripts necessary to reproduce Figures 2, 3, and 7 from the EicC Sullivan process paper (*arXiv:2512.01720*). The primary goal of this project is to achieve strict, non-manipulated replication of the published distributions by accurately enforcing coordinate geometry, absolute cross-section normalizations, and intricate kinematic boundaries.

## Directory Structure
* `tagged-neutron-DIS/`: Contains the main Monte Carlo generation framework for the Sullivan process ($ep \rightarrow e'n\pi^+$), the primary generator scripts (`TaggedN_DIS.cpp`, `run_sim_pion.cpp`), the generated root datasets, and `plot_fig2_3.C`.
* `tagged-Lambda-DIS/`: Contains generic Pythia DIS background simulation data and `plot_fig7.C`, which overlays signal and background distributions.
* `output_plots/`: Contains the final, high-fidelity `.png` reproductions of the paper's figures.

## Detailed Reproduction Steps

### 1. Generating the Sullivan Signal
- **Generator Framework**: The `TaggedN_DIS.cpp` generator simulates the Sullivan process by sampling a specific hyperrectangle in phase space over $\Delta x_B$, $\Delta Q^2$, and $\Delta x_L$.
- **Kinematic Bounding**: We explicitly restrict the simulation generation box to $0.75 < x_L < 1.0$, exactly as stated in the paper. Attempting to expand generation below this point breaks the uniform phase-space normalization because the explicit requirement of $|t| < 1.0 \text{ GeV}^2$ makes generating events below $x_L = 0.363$ mathematically impossible, triggering a hidden rejection loop. 
- **Execution**: The signal is executed via `run_sim_pion.cpp`, which produces $1,000,000$ accepted events saved into `TaggedNeutron-DIS-EicC.root`.

### 2. Physical Cross-Section Normalization
To prevent any manual scaling, arbitrary cutoffs, or "manipulated" visual data, all event weights are physically tied to an integrated luminosity of $1 \text{ pb}^{-1}$.
- **Phase Space Volume**: $V_{sig} = \Delta t \times \Delta Q^2 \times \Delta x_L \times \Delta x_B = 1.0 \times 49.0 \times 0.249 \times 0.99 \approx 12.07 \text{ GeV}^2$.
- **Event Weight Calculation**: Each event in the TTree carries a differential cross section `d4sigma`. The final plotted weight is $w = \frac{L \cdot V_{sig}}{N_{gen}} \cdot d^4\sigma$.
- By strictly adhering to this uniform integration, the Sullivan peak in Figure 7 naturally bounds itself around $\sim 10,000 - 13,000$ events.

### 3. Figure 2: Phase Space Heatmaps
- **Script**: `plot_fig2_3.C`
- **Heatmap Resolution**: The 2D event rate histograms ($x_B$ vs $Q^2$, and $x_\pi$ vs $|t|$) are finely binned at $200 \times 200$. This specific resolution smoothly spreads the rate density natively under $1.0 \text{ Hz}$ per bin, perfectly matching the bounds of the paper's $[10^{-10}, 1]$ Z-axis color scale.
- **Physical Exclusions**: The script explicitly blanks out the region around $x_B \in [0.25, 0.3]$ at lower $Q^2$ by enforcing strict condition cuts: $W^2 > 4 \text{ GeV}^2$, $0.01 < |t| < 1.0 \text{ GeV}^2$, and $x_L > 0.75$.

### 4. Figure 3: Absolute Lab Frame Coordinate Alignment
- **Script**: `plot_fig2_3.C`
- **Crossing Angle Symmetry**: The incident proton beam in the EicC lab frame setup possesses a $50 \text{ mrad}$ ($\approx 2.86^\circ$) crossing angle directed roughly along the $-z$ axis.
- **Neutron Deflection Fix**: Because Pythia standard tracks absolute angles relative to the $+z$ axis ($\sim 177^\circ$), the outgoing neutron polar angle ($\theta_n$) must be reflected via `180.0 - theta_lab`. This precisely translates the neutron distribution so that it symmetrically peaks directly over the $2.86^\circ$ reference line in Figure 3(b).

### 5. Figure 7: Relative Coordinate Geometry
- **Script**: `plot_fig7.C`
- **Beam-Relative Vectors**: Unlike Figure 3, the pseudorapidity ($\eta$) and transverse momentum ($p_T$) in Figure 7 are defined *relative* to the incident proton beam vector. We construct a 3-vector `pBeam_v3` encoding the $50 \text{ mrad}$ angle and use `.Angle(pBeam_v3)` to derive relative scattering angles.
- **Strict Filtering**:
  - **Figure 7(a)**: Cuts enforce $x_L > 0.75$, $M_X > 0.5 \text{ GeV}$, $W^2 > 4 \text{ GeV}^2$, alongside the background baseline $0.01 < |t| < 1.0 \text{ GeV}^2$.
  - **Figure 7(b)**: Cuts enforce $\eta > 5.0$, $p_T < 0.3 \text{ GeV}$, $M_X > 0.5 \text{ GeV}$, $W^2 > 4 \text{ GeV}^2$. The hardcoded Y-axis constraints were removed so the entire un-cut peak shape is visible.
