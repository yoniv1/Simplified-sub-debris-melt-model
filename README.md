# Debris-covered Glacier Surface Mass Balance Model

A physically based surface mass balance (SMB) model for debris-covered glaciers, combining simplified and full surface energy balance (SEB) formulations with explicit debris heat conduction and internal heat storage (IHS).

The model was developed and applied for the Djankuat Glacier (Caucasus) to simulate melt beneath supraglacial debris and investigate the influence of model complexity on sub-debris melt rates and Østrem curves.

---

# Features

- Clean-ice and debris-covered SMB modelling
- Simplified and full surface energy balance options
- Optional internal heat storage within the debris layer
- Optional patchy/thin debris melt enhancement parameterization

---

# Model configurations

The code supports three principal configurations:

| Configuration | Surface energy balance | Internal heat storage |
|---|---|---|
| Full SEB + full IHS | Explicit SEB | Explicit heat conduction |
| Simple SEB + full IHS | Parameterized SEB | Explicit heat conduction |
| Simple SEB + simple IHS | Parameterized SEB | Linear temperature gradient |

# Main model switches

| Switch | Description |
|---|---|
| `full_seb_on` | Enables full surface energy balance |
| `lin_temp_grad_on` | Uses simplified linear temperature gradients |
| `patchy_deb_on` | Enables patchy debris melt enhancement |

---

# Input data

The model requires:

## Meteorological forcing

- Air temperature
- Precipitation
- Incoming shortwave radiation
- Relative humidity (full SEB)
- Wind speed (full SEB)
- Atmospheric emissivity or longwave radiation (full SEB)
- Air pressure (full SEB)

## Glacier data

- DEM
- Glacier mask
- Debris thickness distribution

---

# Citation

If you use this model, please cite:

> Verhaegen et al. (in preparation / submitted)
