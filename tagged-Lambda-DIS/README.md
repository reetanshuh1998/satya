# tagged-Lambda-DIS

Event generator for leading Lambda baryon tagged DIS process ($ep \to e'\Lambda K^+ X$), useful for constraining the **kaon structure function** at EicC.

## Event Generation

**Using ROOT to execute the codes:**
```bash
root -l -b -q 'run_sim_kaon.cpp(1000000)'
# → Creates: TaggedLambda-DIS-EicC.root
```

**Using ACLiC to compile and execute:**
```bash
root -l
[0] .x run_sim_kaon.cpp+
```

## Plotting Scripts

### Figure 7 — Signal vs DIS Background
```bash
root -l -b -q plot_fig7.C
# → Creates: Figure7.png, Figure7_pT03.png, Figure7_pT05.png
```

> **Note:** This script reads the Sullivan signal from `../tagged-neutron-DIS/TaggedNeutron-DIS-EicC-Fig7.root` (generated with $x_L > 0.36$) and the DIS background from `../DIS_Background_EicC.root`. See the main README for why artificial scaling is applied.

### Structure Functions — $F_2^\pi$, $F_2^K$, $xu_v$ vs Drell-Yan Data
```bash
root -l -b -q plot_F2_pion_kaon.C
# → Creates: F2_pion_kaon.png, xuv_pion_vs_E615.png,
#            F2_pion_vs_kaon.png, F2_ratio_kaon_pion.png
```

## PDF Grid Files

| File | Selected via | Content |
|:---|:---:|:---|
| `grid_pion_global_fit.dat` | `setDataSet(211)` | Pion PDFs from global fit (default) |
| `grid_pion_SetA.dat` | `setDataSet(1)` | Pion PDFs — Set A (harder) |
| `grid_pion_SetB.dat` | `setDataSet(2)` | Pion PDFs — Set B |
| `grid_kaon_global_fit.dat` | `setDataSet(321)` | Kaon PDFs from global fit |

## Reference

G. XIE et al., "Tackling the kaon structure function at EicC",
Chinese Physics C 46 (2022) 6, 064107 [arXiv:2109.08483]

Download ROOT: https://root.cern.ch/
Contact: rwang@impcas.ac.cn
