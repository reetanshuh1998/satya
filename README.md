# EicC Sullivan Process Reproduction

## Objective
This repository contains the simulation framework, physics validation, and plotting scripts necessary to reproduce Figures 2, 3, and 7 from the EicC Sullivan process paper (*arXiv:2512.01720*). The primary goal of this project is to achieve strict, non-manipulated replication of the published distributions by accurately enforcing coordinate geometry, absolute cross-section normalizations, and the paper's kinematic selection criteria.

## Directory Structure
* `tagged-neutron-DIS/`: Contains the main Monte Carlo generation framework for the Sullivan process ($ep \rightarrow e'n\pi^+$), the primary generator scripts (`TaggedN_DIS.cpp`, `run_sim_pion.cpp`), the generated root datasets, and `plot_fig2_3.C`.
* `tagged-Lambda-DIS/`: Contains generic Pythia DIS background simulation data and `plot_fig7.C`, which overlays signal and background distributions.
* `output_plots/`: Contains the final, high-fidelity `.png` reproductions of the paper's figures.

## Detailed Reproduction Steps

### 1. Generating the Sullivan Signal
- **Generator Framework**: The `TaggedN_DIS.cpp` generator simulates the Sullivan process by sampling a hyperrectangle in phase space over $\Delta x_B$, $\Delta Q^2$, $\Delta x_L$, and $\Delta t$.
- **Kinematic Range**: The generator samples the **full physical range** $0.001 < x_L < 0.999$ (covering $x_L \in [0, 1]$). Analysis-level cuts such as $x_L > 0.75$ are applied at the plotting stage, not during generation. Generating over the full range is correct because:
  - The Sullivan cross section is weighted per-event by `d4sigma`, so low-$x_L$ events naturally receive negligible weight (the pion flux $f_{\pi/p}(x_L, t)$ peaks near $x_L \approx 1$).
  - It allows Figure 3 to show the complete lab-frame kinematic coverage without pre-selection bias.
  - It mirrors standard physics MC practice: generate inclusively, then apply analysis cuts.
- **Pion Flux Factor**: The pion-exchange flux uses the Bishari form with an exponential regulator: $f_{\pi/p}(x_L, t) = \frac{1}{2\pi} \frac{g^2_{\pi NN}}{4\pi} (1-x_L) \frac{(-t)}{(m_\pi^2 - t)^2} \exp\!\left(-\frac{\Lambda^2(m_\pi^2 - t)}{1 - x_L}\right)$, with $g^2_{\pi NN}/4\pi = 13.6$ and $\Lambda = 0.93~\text{GeV}$. These are standard values from the Sullivan process literature.
- **Execution**: The signal is executed via `run_sim_pion.cpp`, which produces $1{,}000{,}000$ accepted events saved into `TaggedNeutron-DIS-EicC.root`.

### 2. Physical Cross-Section Normalization
All event weights are physically tied to an assumed luminosity; no manual rescaling is applied.
- **Phase Space Volume**: $V = \Delta|t| \times \Delta Q^2 \times \Delta x_L \times \Delta x_B = 0.99 \times 49.0 \times 0.998 \times 0.99 \approx 48.40~\text{GeV}^4$, reflecting the full $x_L \in [0.001, 0.999]$ generation range.
- **Event Weight**: Each event carries a differential cross section `d4sigma`. The plotted weight is $w = \frac{L \cdot V}{N_{\rm gen}} \cdot d^4\sigma$, where $L$ is the integrated luminosity.
- **Figure 7 Luminosity**: $L = 1~\text{pb}^{-1} = 1000~\text{nb}^{-1}$, consistent with the paper's stated running conditions. The Sullivan peak in Figure 7 naturally reaches $\sim 10{,}000$–$13{,}000$ events under this normalization.

### 3. Figure 2: Phase Space Heatmaps
- **Script**: `plot_fig2_3.C`
- **Analysis Cuts**: The paper-defined selection $x_L > 0.75$, $0.01 < |t| < 1.0~\text{GeV}^2$, $W^2 > 4~\text{GeV}^2$, $1 < Q^2 < 50~\text{GeV}^2$ is enforced at fill time, not at generation time.

### 4. Figure 3: Absolute Lab Frame Coordinate Alignment
- **Script**: `plot_fig2_3.C`
- **Neutron Angle**: The neutron polar angle $\theta_n$ is measured from the $-z$ axis (the nominal ion beam direction) via `180° - theta_ROOT`, where `theta_ROOT` is ROOT's angle from $+z$. This correctly places the distribution peak at the $50~\text{mrad}~(\approx 2.86^\circ)$ beam crossing angle.

### 5. Figure 7: Signal/Background Separation and pT Comparison
- **Script**: `plot_fig7.C`
- **Beam-Relative Kinematics**: The neutron $\eta$ and $p_T$ in Figure 7 are computed relative to the incident proton beam direction (50 mrad crossing angle) using `.Angle(pBeam_v3)`.
- **Two pT Cases**: The paper has an internal inconsistency — the body text and figure caption state $p_T^n < 0.3~\text{GeV}$, while the label embedded in the figure reads "$p_T < 0.5~\text{GeV}$". The script therefore produces **two outputs**:
  - `Figure7_pT03.png`: $p_T^n < 0.3~\text{GeV}$ (matches caption/body text).
  - `Figure7_pT05.png`: $p_T^n < 0.5~\text{GeV}$ (matches the in-figure label).
