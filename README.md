# REE cross-plot: MREE/MREE\* vs HREE/LREE

An interactive R Shiny tool that places your rare earth element (REE) data on the MREE/MREE\* vs HREE/LREE cross-plot (idea after Du et al., 2016; Blaser et al., 2020, Fig. 3), together with

- the marine **authigenic array** (seawater to Fe-Mn nodule mixing),
- the **detrital array** (MORB, Greenland Archaean crust and Canadian Precambrian shield mixing),
- kernel-density **fields** of published foraminifera, pore water and seawater data (including GEOTRACES IDP2025 seawater, split at 200 m).

The app runs entirely in the browser (Shinylive / WebAssembly), so uploaded data never leave your computer.

## Using the app

The sidebar has six collapsible sections.

1. **Your data**: upload a CSV/TSV (or tick *Use example data*; uploading a file unticks it), choose concentrations or ratios, the unit basis (mass or molar), ID and group columns, and further columns to carry into the results. The columns recognised as REE (and Y) are listed under the status message.
2. **Normalisation & ratios**: normalising reference (PAAS Pourmand et al. 2012 by default; PAAS McLennan 1989, NASC, EUS Bau et al. 2018, CI chondrites after Pourmand 2012 or O'Neill 2016, or your own one-row reference file), LREE definition (La-Pr-Nd or Pr-Nd-Sm) gap method (interpolate, group mean, none) and the Ce/Ce\* definition.
3. **Screening & symbols**: symbol fill by group, a numeric column (e.g. age, depth), REE roughness or Gd/Gd\*; an optional roughness limit; a minimum number of measured REE and required elements; group visibility (tick a single group to show only that one, or untick groups to exclude them). Screening applies to your data only.
4. **Plot**: λ space (experimental), log10 axes (both axes together), legend on/off, plot window, field coverage and smoothness, point size and optional density fields for your own groups.
5. **Reference layers**: density fields, reference points, arrays, labels, and the individual reference datasets.
6. **Anthropogenic Gd**: use or ignore measured Gd, the method for natural Gd (Bau and Dulski 1996, its geometric version, Sm²/Nd, Tb²/Dy, or BLambdaR, experimental) and the natural Gd/Gd\* baseline.

Options that have no effect in λ space are greyed out and annotated while λ space is ticked.

The **Spectra** tab plots the REE patterns of your data (normalised to the chosen reference, optionally to one element; element or ionic-radius axis after Shannon, 1976; Y optional), with the same samples and group colours as the cross-plot.

Downloads: figure (PDF; PNG where supported), results table (CSV), arrays and fields in the current configuration (CSV), spectra (PDF) and the normalised REE (CSV). The results table lists every row of the upload in the original order; the column `row` is the row number in your file (header not counted) and identifies a sample uniquely. The calculation settings are written into the table. Click a point to see its values.

## Methods in brief

| Quantity | Definition (shale-normalised values) |
|---|---|
| HREE/LREE | (Tm + Yb + Lu) / LREE, with LREE = La + Pr + Nd (default) or Pr + Nd + Sm |
| MREE/MREE\* | (Gd + Tb + Dy) / ((HREE + LREE) / 2) |
| Anomalies (table only) | geometric, after Barrat and Bayon (2024): Ce/Ce\* = Ce·Nd / Pr² (arithmetic 2 Ce / (La + Pr) and Ce / sqrt(La Pr) selectable), La/La\* = La·Nd² / Pr³, Eu/Eu\* = Eu / sqrt(Sm·Gd) |
| Gd/Gd\* | Gd / Gd\*, with Gd\* = 0.33 Sm + 0.67 Tb (Bau and Dulski, 1996) by default; alternatives Sm^(1/3) Tb^(2/3), Sm²/Nd, Tb²/Dy (Barrat and Bayon, 2024) or BLambdaR |
| Anthropogenic Gd | Gd · max(0, 1 − baseline / (Gd/Gd\*)), in the input units |
| REE roughness | square root of the summed squared deviations of Nd, Tb, Dy, Ho, Er, Tm and Yb from the mean of their neighbours (unpublished screening measure; needs all 11 elements measured) |
| λ space | −λ1 and −λ2 from a base-R port of BLambdaR (Anenburg and Williams, 2022; O'Neill, 2016 radii); computed only when ticked |

**Units.** Mass data are normalised to mass-based standards, molar data (pmol/kg) to the standard divided by atomic mass. Mixing the two biases HREE/LREE by about 18%.

**Gaps.** Missing Pr, Tb, Ho, Er, Tm, Lu are placed on the log-linear line through their neighbours as a function of atomic number (Pm omitted, Barrat and Bayon, 2024), or replaced by group means, or left empty.

**Arrays and fields.** Arrays are linear mixtures in normalised space, recomputed for the chosen standard and LREE definition. Fields are 2-D kernel density estimates (`MASS::kde2d`) enclosing the chosen share of the data. Reference layers always use measured Gd and interpolated gaps; only the standard and LREE definition follow your choice.

To reproduce Fig. 3 of Blaser et al. (2020), choose PAAS McLennan (1989) and untick the GEOTRACES datasets.

## Verification

`tests/run_tests.R` (91 checks) covers the core calculations, lambda module, arrays, fields and parsing; `tests/test_app.R` (41 checks) exercises the Shiny server logic. From the repository root:

```sh
export LC_ALL=C.UTF-8
Rscript tests/run_tests.R
Rscript tests/test_app.R
Rscript tests/render_check.R figure.png   # args: out.png [norm] [lree] [space] [log] [example]
```

## Run locally

```r
install.packages(c("shiny", "ggplot2", "MASS"))
shiny::runApp("app")
```

## Deploy on GitHub Pages and embed on a website

1. Create an empty public repository `ree-crossplot` on GitHub, then from this folder: `git remote add origin git@github.com:<user>/ree-crossplot.git` and `git push -u origin main`.
2. In the repository, *Settings > Pages > Source*: **GitHub Actions**.
3. The workflow `.github/workflows/deploy.yml` runs the tests, exports the app with `shinylive::export("app", "_site")` and publishes it at `https://<user>.github.io/ree-crossplot/` (check the *Actions* tab; start it by hand with *Run workflow* if needed). Open the URL and test the app once; the first load downloads the R environment (WebAssembly) and takes a while.
4. In the repository of the website (Beautiful Jekyll), copy `embed/REE-crossplot.md` to `Other/REE-crossplot.md`, replace `YOUR-GITHUB-USERNAME` (three places) and add the page to the navigation in `_config.yml`:

```yaml
  Other:
   - Incomplete Paleo Data Demo: "Other/PaleoDataDemo"
   - REE cross-plot: "Other/REE-crossplot"
```

5. Commit and push the website; Jekyll rebuilds within a minute or two. Later changes to the app only need a push to the app repository.

## Rebuilding the reference data

`Rscript data-raw/build_data.R` rebuilds `app/data/` from the compilation files and the GEOTRACES subset in `data-raw/` (GEOTRACES: quality flag "no problem reported", La and depth present; 118 rows duplicating compilation seawater removed).

## Limitations

- The Shinylive build has not been tested in a browser (only R code, server logic and the native Shiny app were tested); run the workflow once and check the page. The fields are recomputed on each change and may feel slower in WebAssembly.
- PNG download may be unavailable in webR; PDF is the safe option.
- The lambda space, λ roughness (default limit 0.03 is provisional), REE roughness and the BLambdaR Gd option are experimental.
- Upload size follows the Shiny default (5 MB); select the REE columns before uploading large files.
- Y of the shale standards (spectra tab only) is taken from the author's REE_spectrum() function and not independently verified; molar data are taken as pmol/kg on the spectra tab.
- The plot is a diagnostic, not a source apportionment.

## Licence

MIT, see `LICENSE`.

## References

- Abbott, A.N., Haley, B.A., McManus, J., Reimers, C.E., 2015. The sedimentary flux of dissolved rare earth elements to the ocean. Geochim. Cosmochim. Acta 154, 186–200. https://doi.org/10.1016/j.gca.2015.01.010
- Anenburg, M., Williams, M.J., 2022. Quantifying the tetrad effect, shape components, and Ce–Eu–Gd anomalies in rare earth element patterns. Math. Geosci. 54, 47–70. https://doi.org/10.1007/s11004-021-09959-5
- Barrat, J.-A., Bayon, G., Tripier, R., Rouget, M.-L., Ben Salem, D., 2026. Anthropogenic gadolinium contamination in tap water across the Paris megacity: Rare earth elements trace sources and mixing. Environ. Pollut. 405, 128606. https://doi.org/10.1016/j.envpol.2026.128606
- Barrat, J.-A., Bayon, G., 2024. Practical guidelines for representing and interpreting rare earth abundances in environmental and biological studies. Chemosphere 352, 141487. https://doi.org/10.1016/j.chemosphere.2024.141487
- Bau, M., Dulski, P., 1996. Anthropogenic origin of positive gadolinium anomalies in river waters. Earth Planet. Sci. Lett. 143, 245–255. https://doi.org/10.1016/S0012-821X(96)00127-6
- Bau, M., Schmidt, K., Koschinsky, A., Hein, J., Kuhn, T., Usui, A., 2014. Discriminating between different genetic types of marine ferro-manganese crusts and nodules based on rare earth elements and yttrium. Chem. Geol. 381, 1–9. https://doi.org/10.1016/j.chemgeo.2014.05.004
- Bau, M., Schmidt, K., Pack, A., Bendel, V., Kraemer, D., 2018. The European Shale: An improved data set for normalisation of rare earth element and yttrium concentrations in environmental and biological samples from Europe. Appl. Geochem. 90, 142–149. https://doi.org/10.1016/j.apgeochem.2018.01.008
- Blaser, P., Gutjahr, M., Pöppelmeier, F., Frank, M., Kaboth-Bahr, S., Lippold, J., 2020. Labrador Sea bottom water provenance and REE exchange during the past 35,000 years. Earth Planet. Sci. Lett. 542, 116299. https://doi.org/10.1016/j.epsl.2020.116299
- Crocket, K.C., Hill, E., Abell, R.E., Johnson, C., Gary, S.F., Brand, T., Hathorne, E.C., 2018. Rare Earth Element Distribution in the NE Atlantic: Evidence for Benthic Sources, Longevity of the Seawater Signal, and Biogeochemical Cycling. Front. Mar. Sci. 5, 147. https://doi.org/10.3389/fmars.2018.00147
- Du, J., Haley, B.A., Mix, A.C., 2016. Neodymium isotopes in authigenic phases, bottom waters and detrital sediments in the Gulf of Alaska and their implications for paleo-circulation reconstruction. Geochim. Cosmochim. Acta 193, 14–35. https://doi.org/10.1016/j.gca.2016.08.005
- Filippova, A., Frank, M., Kienast, M., Rickli, J., Hathorne, E., Yashayaev, I.M., Pahnke, K., 2017. Water mass circulation and weathering inputs in the Labrador Sea based on coupled Hf–Nd isotope compositions and rare earth element distributions. Geochim. Cosmochim. Acta 199, 164–184. https://doi.org/10.1016/j.gca.2016.11.024
- Gale, A., Dalton, C.A., Langmuir, C.H., Su, Y., Schilling, J.-G., 2013. The mean composition of ocean ridge basalts. Geochem. Geophys. Geosyst. 14, 489–518. https://doi.org/10.1029/2012GC004334
- GEOTRACES Intermediate Data Product Group, 2025. The GEOTRACES Intermediate Data Product 2025 (IDP2025). https://doi.org/10.5285/42C92148-8D03-8BE6-E063-7086ABC09F0C
- Hathorne, E.C., Stichel, T., Brück, B., Frank, M., 2015. Rare earth element distribution in the Atlantic sector of the Southern Ocean: The balance between particle scavenging and vertical supply. Mar. Chem. 177, 157–171. https://doi.org/10.1016/j.marchem.2015.03.011
- Kulaksız, S., Bau, M., 2007. Contrasting behaviour of anthropogenic gadolinium and natural rare earth elements in estuaries and the gadolinium input into the North Sea. Earth Planet. Sci. Lett. 260, 361–371. https://doi.org/10.1016/j.epsl.2007.06.016
- Martin, E.E., Blair, S.W., Kamenov, G.D., Scher, H.D., Bourbon, E., Basak, C., Newkirk, D.N., 2010. Extraction of Nd isotopes from bulk deep sea sediments for paleoceanographic studies on Cenozoic time scales. Chem. Geol. 269, 414–431. https://doi.org/10.1016/j.chemgeo.2009.10.016
- McLennan, S.M., 1989. Chapter 7. Rare earth elements in sedimentary rocks: influence of provenance and sedimentary processes, in: Geochemistry and Mineralogy of Rare Earth Elements. De Gruyter, Berlin, Boston. https://doi.org/10.1515/9781501509032-010
- Nance, W.B., Taylor, S.R., 1976. Rare earth element patterns and crustal evolution—I. Australian post-Archean sedimentary rocks. Geochim. Cosmochim. Acta 40, 1539–1551.
- O'Neill, H.St.C., 2016. The smoothness and shapes of chondrite-normalized rare earth element patterns in basalts. J. Petrol. 57, 1463–1508. https://doi.org/10.1093/petrology/egw047
- Osborne, A.H., Hathorne, E.C., Schijf, J., Plancherel, Y., Böning, P., Frank, M., 2017. The potential of sedimentary foraminiferal rare earth element patterns to trace water masses in the past. Geochem. Geophys. Geosyst. 18, 1550–1568. https://doi.org/10.1002/2016GC006782
- Palmer, M.R., 1985. Rare earth elements in foraminifera tests. Earth Planet. Sci. Lett. 73, 285–298. https://doi.org/10.1016/0012-821X(85)90077-9
- Pourmand, A., Dauphas, N., Ireland, T.J., 2012. A novel extraction chromatography and MC-ICP-MS technique for rapid analysis of REE, Sc and Y: Revising CI-chondrite and Post-Archean Australian Shale (PAAS) abundances. Chem. Geol. 291, 38–54. https://doi.org/10.1016/j.chemgeo.2011.08.011
- Shannon, R.D., 1976. Revised effective ionic radii and systematic studies of interatomic distances in halides and chalcogenides. Acta Crystallogr. A 32, 751–767. https://doi.org/10.1107/S0567739476001551
- Shaw, D.M., Cramer, J.J., Higgins, M.D., Truscott, M.G., 1986. Composition of the Canadian Precambrian shield and the continental crust of the earth. Geol. Soc. Lond. Spec. Publ. 24, 275–282. https://doi.org/10.1144/GSL.SP.1986.024.01.24
- Taylor, S.R., McLennan, S.M., 1985. The Continental Crust: Its Composition and Evolution. Blackwell, Oxford.
- Wedepohl, K.H., Heinrichs, H., Bridgwater, D., 1991. Chemical characteristics and genesis of the quartz-feldspathic rocks in the Archean crust of Greenland. Contrib. Mineral. Petrol. 107, 163–179. https://doi.org/10.1007/BF00310705
- Zheng, X.-Y., Plancherel, Y., Saito, M.A., Scott, P.M., Henderson, G.M., 2016. Rare earth elements (REEs) in the tropical South Atlantic and quantitative deconvolution of their non-conservative behavior. Geochim. Cosmochim. Acta 177, 217–237. https://doi.org/10.1016/j.gca.2016.01.018
