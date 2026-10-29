---
layout: page
title: REE cross-plot (MREE/MREE* vs HREE/LREE)
subtitle: Interactive tool for rare earth element patterns in paleoceanography
---

This tool places your rare earth element (REE) data on the MREE/MREE\* vs HREE/LREE cross-plot, together with the marine authigenic array, the detrital array and density fields of published foraminifera, pore water and seawater data (including GEOTRACES IDP2025). The diagram follows Du et al. (2016) and Blaser et al. (2020, Fig. 3).

**How to use it.** Tick *Use example data* or upload a CSV/TSV with one row per sample and one column per element (`La`, `Ce`, ..., `Lu`; names such as `La_ppm` work). Choose mass or molar units, the shale or chondrite normalisation, and optionally a group column. The *Spectra* tab shows the normalised REE patterns, *Computed values* lists all results, and the buttons below the plot download the figure and the data. Your data are processed in your browser and are not uploaded to any server.

The first load takes a little while, because the R environment (WebAssembly) is downloaded. For a larger view, [open the tool in its own window](https://YOUR-GITHUB-USERNAME.github.io/ree-crossplot/){:target="_blank"}.

<iframe
  src="https://YOUR-GITHUB-USERNAME.github.io/ree-crossplot/"
  title="REE cross-plot: MREE/MREE* vs HREE/LREE"
  style="width:100%; height:1150px; border:1px solid #ccc; border-radius:6px;"
  loading="lazy"
  allow="clipboard-write"></iframe>

**Source code and documentation:** [github.com/YOUR-GITHUB-USERNAME/ree-crossplot](https://github.com/YOUR-GITHUB-USERNAME/ree-crossplot). The methods, options and all references are in the *About & references* tab of the tool.

**Please cite** Blaser et al. (2020) and Du et al. (2016) when using the diagram:

- Blaser, P., Gutjahr, M., Pöppelmeier, F., Frank, M., Kaboth-Bahr, S., Lippold, J., 2020. Labrador Sea bottom water provenance and REE exchange during the past 35,000 years. Earth Planet. Sci. Lett. 542, 116299. https://doi.org/10.1016/j.epsl.2020.116299
- Du, J., Haley, B.A., Mix, A.C., 2016. Neodymium isotopes in authigenic phases, bottom waters and detrital sediments in the Gulf of Alaska and their implications for paleo-circulation reconstruction. Geochim. Cosmochim. Acta 193, 14–35. https://doi.org/10.1016/j.gca.2016.08.005
