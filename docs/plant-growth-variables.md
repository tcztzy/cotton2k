# Plant growth and root variables (Rust port)

This note documents the intent and equations behind the new variables introduced in the Rust
translations of the plant growth and root growth modules. The equations below follow the
implementation exactly and mirror the original model logic.

## Potential leaf growth (`potential_leaf_growth`)

Key variables
- `V_RATIO`: ratio of actual leaf + petiole growth to their potential requirements (set in
  `dry_matter_balance`, used in `actual_leaf_growth`).
- `LeafWeightAreaRatio`: dry matter per unit leaf area (g/dm^2).
- `PotGroLeafArea*`, `PotGroLeafWeight*`, `PotGroPetioleWeight*`: potential area/weight increments
  for leaves and petioles on prefruiting nodes, main stem nodes, and fruiting nodes.
Notation: `v_i` in this section refers to `VPOTLF[i]`.

Monomolecular leaf area growth
```{math}
L(t) = s_{\max} \left(1 - e^{-c t^{p}}\right)
```
```{math}
\frac{dL}{dt} = r = s_{\max} c p e^{-c t^{p}} t^{p-1}
```

Water stress reduction and leaf weight/area ratio
```{math}
\mathrm{wstrlf} = \mathrm{WaterStress} \left(1 + v_0 (2 - \mathrm{WaterStress})\right) - v_0
```
```{math}
\mathrm{wtfstrs} = v_1 + v_2 (1 - \mathrm{wstrlf})
```
```{math}
\mathrm{LeafWeightAreaRatio} = \frac{\mathrm{wtfstrs}}{v_4 + T (v_5 - T v_6)}
```
Where `T` is the limited daily mean temperature (`AvrgDailyTemp` but not below `v_3`).

Potential area and weight increments (example for a given leaf)
```{math}
\mathrm{PotGroLeafArea} = r \; \mathrm{wstrlf} \; \mathrm{pixda} \; f_T(T)
```
```{math}
\mathrm{PotGroLeafWeight} = \mathrm{PotGroLeafArea} \; \mathrm{LeafWeightAreaRatio}
```
```{math}
\mathrm{PotGroPetioleWeight} = \mathrm{PotGroLeafArea} \; \mathrm{LeafWeightAreaRatio} \; v_{13}
```
`f_T(T)` is the temperature response function defined below.

## Temperature response for leaf growth (`temperature_on_leaf_growth_rate`)

Piecewise quadratic (normalized by `PAR[7]`):
Notation: `PAR_i` refers to `PAR[i]`.
```{math}
ra(T) =
\begin{cases}
PAR_1 + T (PAR_2 - T PAR_3), & T > PAR_0 \\
PAR_4 + T (PAR_5 - T PAR_6), & T \le PAR_0
\end{cases}
```
```{math}
f_T(T) = \max(0, ra(T)) / PAR_7
```

## Dry matter balance and allocation (`dry_matter_balance`)

Notation: `v_i` in this section refers to `VCHBAL[i]`. `N_fr`, `N_veg`, `N_root` map to
`NStressFruiting`, `NStressVeg`, `NStressRoots`.

Carbohydrate demand (g plant^-1 day^-1):
```{math}
cd_{sqar} = \mathrm{PotGroAllSquares} \frac{N_{fr} + v_0}{v_0 + 1}
```
```{math}
cd_{boll} = (\mathrm{PotGroAllBolls} + \mathrm{PotGroAllBurrs}) \frac{N_{fr} + v_0}{v_0 + 1}
```
```{math}
cd_{leaf} = \mathrm{PotGroAllLeaves} \frac{N_{veg} + v_1}{v_1 + 1}
```
```{math}
cd_{stem} = \mathrm{PotGroStem} \frac{N_{veg} + v_2}{v_2 + 1}
```
```{math}
cd_{root} = \mathrm{PotGroAllRoots} \frac{N_{root} + v_3}{v_3 + 1}
```
```{math}
cd_{pet} = \mathrm{PotGroAllPetioles} \frac{N_{veg} + v_{14}}{v_{14} + 1}
```
```{math}
cd_{sum} = cd_{stem} + cd_{leaf} + cd_{pet} + cd_{root} + cd_{sqar} + cd_{boll}
```

Available carbohydrates and carbon stress:
```{math}
C_{pool} = \mathrm{NetPhotosynthesis} + \mathrm{ReserveC} \; v_{13}
```
```{math}
\mathrm{CarbonStress} = \min\left(1, \frac{C_{pool}}{cd_{sum}}\right)
```

Allocation under limitation (`CarbonStress < 1`):
```{math}
bsratio = \frac{C_{pool}}{cd_{boll} + cd_{sqar}}
```
```{math}
ffr = (v_5 + v_6 (1 - \mathrm{WaterStress})) \; bsratio
```
```{math}
ffr \leftarrow \mathrm{clamp}(ffr, 0, \min(1, bsratio))
```
```{math}
pd_{boll} = cd_{boll} \; ffr, \quad pd_{sqar} = cd_{sqar} \; ffr
```
```{math}
C_{avail} = C_{pool} - pd_{boll} - pd_{sqar}
```
```{math}
flf = v_7 \frac{C_{avail}}{cd_{leaf} + cd_{pet}}, \quad flf \leftarrow \mathrm{clamp}(flf, 0, 1)
```
```{math}
\mathrm{TotalActualLeafGrowth} = cd_{leaf} \; flf
```
```{math}
\mathrm{TotalActualPetioleGrowth} = cd_{pet} \; flf
```
```{math}
C_{avail} \leftarrow C_{avail} - \mathrm{TotalActualLeafGrowth} - \mathrm{TotalActualPetioleGrowth}
```
```{math}
ratio = v_8 + v_9 e^{-v_{10} (W_{stem} + W_{leaf} + W_{pet}) \; \mathrm{PerPlantArea}}
```
```{math}
ratio \leftarrow ratio \; v_{11}
```
```{math}
rt_{max} = \frac{ratio}{ratio + 1} \left(1 + v_{12} (1 - \mathrm{WaterStress})\right)
```
```{math}
frt = rt_{max} \frac{C_{avail}}{cd_{root}}, \quad frt \leftarrow \mathrm{clamp}(frt, 0, 1)
```
```{math}
\mathrm{CarbonAllocatedForRootGrowth} = \max(cd_{root} \; frt, C_{avail} - cd_{stem})
```
```{math}
C_{avail} \leftarrow C_{avail} - \mathrm{CarbonAllocatedForRootGrowth}
```
```{math}
fst = \frac{C_{avail}}{cd_{stem}}, \quad fst \leftarrow \mathrm{clamp}(fst, 0, 1)
```
```{math}
\mathrm{ActualStemGrowth} = cd_{stem} \; fst
```
```{math}
\mathrm{ExtraCarbon}_1 = \max(0, C_{avail} - \mathrm{ActualStemGrowth})
```

Reserve update and derived ratios:
```{math}
\mathrm{ReserveC} \leftarrow \mathrm{ReserveC} + \mathrm{NetPhotosynthesis} -
(\mathrm{ActualStemGrowth} + \mathrm{TotalActualLeafGrowth} + \mathrm{TotalActualPetioleGrowth} +
\mathrm{CarbonAllocatedForRootGrowth} + pd_{boll} + pd_{sqar})
```
```{math}
\mathrm{ResMax} = v_4 \; \mathrm{TotalLeafWeight}
```
```{math}
\mathrm{ExtraCarbon}_2 = \max(0, \mathrm{ReserveC} - \mathrm{ResMax})
```
```{math}
\mathrm{ReserveC} \leftarrow \min(\mathrm{ReserveC}, \mathrm{ResMax})
```
```{math}
\mathrm{ExtraCarbon} = \mathrm{ExtraCarbon}_1 + \mathrm{ExtraCarbon}_2
```
```{math}
\mathrm{FruitGrowthRatio} = \frac{pd_{sqar} + pd_{boll}}{\mathrm{PotGroAllSquares} + \mathrm{PotGroAllBolls} + \mathrm{PotGroAllBurrs}}
```
```{math}
V\_RATIO = \frac{\mathrm{TotalActualLeafGrowth} + \mathrm{TotalActualPetioleGrowth}}{\mathrm{PotGroAllLeaves} + \mathrm{PotGroAllPetioles}}
```

## Actual fruit growth (`actual_fruit_growth`)

Actual growth uses `FruitGrowthRatio` as a proportional scaler:
```{math}
\Delta W_{square} = \mathrm{PotGroSquares} \; \mathrm{FruitGrowthRatio}
```
```{math}
\Delta W_{boll} = \mathrm{PotGroBolls} \; \mathrm{FruitGrowthRatio}
```
```{math}
\Delta W_{burr} = \mathrm{PotGroBurrs} \; \mathrm{FruitGrowthRatio}
```

## Actual leaf growth (`actual_leaf_growth`)

Actual leaf and petiole growth are scaled by `V_RATIO`:
```{math}
\Delta W_{leaf} = \mathrm{PotGroLeafWeight} \; V\_RATIO
```
```{math}
\Delta W_{pet} = \mathrm{PotGroPetioleWeight} \; V\_RATIO
```
```{math}
\Delta A_{leaf} = \mathrm{PotGroLeafArea} \; V\_RATIO
```

## Dry matter balance check (`check_dry_matter_balance`)

```{math}
\mathrm{avail} = \mathrm{PlantWeightAtStart} + \mathrm{CumNetPhotosynth}
```
```{math}
\mathrm{PlantWeight} = W_{root} + W_{stem} + W_{leaf} + W_{pet} + W_{square} +
W_{boll,green} + W_{boll,open} + W_{burr,green} + W_{burr,open} + \mathrm{ReserveC}
```
```{math}
\mathrm{used} = \mathrm{PlantWeight} + \mathrm{GreenBollsLost} + \mathrm{AbscisedLeafWeight} +
\mathrm{BloomWeightLoss} + \mathrm{RootWeightLoss}
```
```{math}
\mathrm{chobal} = \mathrm{avail} - \mathrm{used}
```

## Defoliation (`defoliate`)

Persistent state variables
- `DEFKGH`: defoliant intercepted by the plant canopy (kg/ha).
- `TDFKGH`: cumulative defoliant amount (kg/ha).
- `IDSW`: switch indicating whether the predicted defoliation date has been set.
Notation: `p_i` corresponds to `P1..P7` constants in the implementation.

When a defoliation is predicted (`DefoliantAppRate[i] < -99`):
```{math}
TDFKGH = 2.5
```

Otherwise, intercepted defoliant is accumulated:
```{math}
DEFKGH \leftarrow DEFKGH + \mathrm{DefoliantAppRate}_i \cdot 0.95 \cdot 1.12085 \cdot 0.75
```
```{math}
TDFKGH \leftarrow TDFKGH + DEFKGH
```
If `DefoliationMethod[i] != 0`, replace `0.95` by `LightIntercept`.

Percent defoliation (bounded to [0, 40]):
```{math}
\mathrm{dum} = -\mathrm{LwpMin} \cdot 10
```
```{math}
\mathrm{PercentDefoliation} = p_1 + p_2 T + p_3 TDFKGH + p_4 (D - D_{first}) +
 p_5 \mathrm{dum} - p_6 \mathrm{dum}^2 + p_7 T \; TDFKGH \; (D - D_{first}) \; \mathrm{dum}
```
where `T = AvrgDailyTemp`, `D = Daynum`.

## Root growth modifiers

Notation: `p_i` refers to each function's local constants; `psi` is soil water potential,
`theta` is volumetric water content.

Soil mechanical resistance (`soil_mechanic_resistance`):
```{math}
rtimpd_{min} = \min\big(\mathrm{RootImpede}_{l,k},\; \mathrm{RootImpede}_{l,k-1},\; \mathrm{RootImpede}_{l,k+1},\; \mathrm{RootImpede}_{l+1,k}\big)
```
```{math}
rtpct = \mathrm{clamp}(p_1 - p_2 \; rtimpd_{min},\; p_3,\; 1)
```

Soil aeration effect (`soil_air_on_root_growth`):
```{math}
rt_{rdo} =
\begin{cases}
 p_2, & \psi > p_1 \\
 1, & \text{otherwise}
\end{cases}
```
```{math}
rt_{rdo} = p_3 \quad \text{if } \theta \ge \mathrm{PoreSpace}
```

Soil nitrate effect (`soil_nitrate_on_root_growth`):
```{math}
rt_{rdn} =
\begin{cases}
 p_2, & \mathrm{VolNo3NContent} < p_1 \\
 1, & \text{otherwise}
\end{cases}
```

Soil water effect (`soil_water_on_root_growth`):
```{math}
smf = \mathrm{clamp}\left(\left(\frac{p_1 + \psi}{p_2}\right)^3,\; 0.02,\; 1\right)
```

## Lateral root initiation and root summation

Lateral root initiation (`initiate_lateral_roots`):
```{math}
sdl \leftarrow \mathrm{TapRootLength} - \mathrm{DepthLastRootLayer}
```
```{math}
\text{for } l = \mathrm{LastTaprootLayer} \ldots 0: \; sdl \leftarrow sdl + dl_l
```
```{math}
\text{if } sdl > \mathrm{DISTLR} \text{ and } \mathrm{LateralRootFlag}_l = 1,\; \mathrm{LateralRootFlag}_l \leftarrow 2
```

Root summation (`root_summation`):
```{math}
roots = \sum_{l,k} \sum_{i=0}^{2} \mathrm{RootWeight}_{l,k,i}
```
```{math}
\mathrm{TotalRootWeight} = roots \cdot 100 \cdot \frac{\mathrm{PerPlantArea}}{\mathrm{RowSpace}}
```
