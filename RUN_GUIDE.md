# Reproduction Guide: EicC Sullivan Process

This guide provides step-by-step instructions on how to generate the Monte Carlo data and produce the exact plots (Figures 2, 3, and 7) for the paper-faithful reproduction of the EicC Sullivan process (arXiv:2512.01720).

> [!NOTE]
> All steps assume you are starting from the root directory of the repository (`/home/reet/satyajit`).

## Step 1: Generate the Signal Data (Sullivan Process)

To correctly model the phase space for all figures, we use two separate signal generations. The paper specifies $x_L > 0.75$ for Figure 2, but Figure 7 requires plotting the discriminator shape down to the physical boundary ($x_L \approx 0.36$).

1. Navigate to the generator directory:
   ```bash
   cd tagged-neutron-DIS
   ```
2. Generate the high-$x_L$ sample for Figures 2 and 3:
   ```bash
   root -l -b -q run_sim_pion.cpp
   ```
   *This saves `TaggedNeutron-DIS-EicC.root`.*
3. Generate the full-$x_L$ sample for Figure 7:
   ```bash
   root -l -b -q run_sim_pion_fig7.cpp
   ```
   *This saves `TaggedNeutron-DIS-EicC-Fig7.root`.*

## Step 2: Generate the Background Data (PYTHIA8 DIS)

The generic deep inelastic scattering (DIS) background is generated using PYTHIA8. The script extracts the exact PYTHIA cross-section (`sigmaGen`) to correctly normalize the background to $1 \text{ pb}^{-1}$.

1. Navigate to the repository root:
   ```bash
   cd /home/reet/satyajit
   ```
2. Compile the generator (you must link both ROOT and PYTHIA8). Depending on your system's PYTHIA8 installation, the include/lib paths might differ. A standard compilation command is:
   ```bash
   g++ -O2 scripts/generate_dis_background.cpp -o generate_dis_background $(root-config --cflags --libs) -I/home/reet/pythia/pythia8312/include -L/home/reet/pythia/pythia8312/lib -lpythia8
   ```
   *(We have updated this command to point directly to your local PYTHIA8 installation at `/home/reet/pythia/pythia8312`)*
3. Run the compiled executable from the repository root:
   ```bash
   export LD_LIBRARY_PATH=/home/reet/pythia/pythia8312/lib:$LD_LIBRARY_PATH
   ./generate_dis_background
   ```
   *This generates 300,000 events and saves the output to `DIS_Background_EicC.root` in the root repository directory.*

> [!IMPORTANT]
> Ensure that `DIS_Background_EicC.root` is located in the root repository directory, as the plotting macros expect it at the `../DIS_Background_EicC.root` relative path.

## Step 3: Reproduce Figures 2 & 3

With the signal data generated, you can now plot the event rate distributions and kinematic scatter plots.

1. Navigate to the signal directory:
   ```bash
   cd tagged-neutron-DIS
   ```
2. Run the plotting macro:
   ```bash
   root -l -b -q plot_fig2_3.C
   ```
3. **Outputs generated:**
   - `Figure2.png` (Event rate in $x_B$-$Q^2$ and $x_\pi$-$|t|$ planes)
   - `Figure3.png` (Electron $\eta$-$p$ and neutron $\theta$-$p$ distributions)

## Step 4: Reproduce Figure 7

Figure 7 overlays the Sullivan process signal against the generic DIS background, demonstrating the suppression power of the geometric and kinematic cuts ($|t| < 1.0$, $x_L > 0.75$, $\eta_n > 5.0$).

1. Navigate to the canonical Figure 7 plotting directory:
   ```bash
   cd ../tagged-Lambda-DIS
   ```
2. Run the plotting macro:
   ```bash
   root -l -b -q plot_fig7.C
   ```
3. **Outputs generated:**
   - `Figure7_pT03.png` / `Figure7.png` (Comparison with $p_T < 0.3$ GeV cut)
   - `Figure7_pT05.png` (Comparison with $p_T < 0.5$ GeV cut)

> [!TIP]
> The plots will display yields corresponding to an integrated luminosity of $1 \text{ pb}^{-1}$, using purely paper-faithful normalization calculations with no empirical scaling factors.
