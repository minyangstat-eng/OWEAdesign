# Introduction slides

A talk that introduces optimal experimental design, the Optimal Weights
Exchange Algorithm (OWEA) of Yang, Biedermann and Tang (2013), the
`OWEAdesign` package and its web app. Use the deck as it is, or adapt it: the
source and every figure are here.

| File | What it is |
|---|---|
| `owea_slides.pdf` | the compiled deck (16:9) |
| `owea_slides.tex` | the LaTeX (beamer, metropolis theme) source |
| `img/` | the figures: screenshots of the app, one per step, and the diagrams |
| `owea_fig.py`, `cmp_fig.py` | Python scripts that draw the diagrams in `img/` |
| `fonts/` | the Lato and DejaVu Sans Mono fonts the deck uses, so it looks the same everywhere |
| `owea_installation.pdf` | a one-page guide to installing R and the package |

## Compiling the deck

```
xelatex owea_slides.tex     # run twice (page numbers, progress bar)
```

`pdflatex` works too (the bundled fonts are then replaced by the default
sans-serif). The metropolis theme is part of TeX Live; nothing else is needed.
After editing a figure script, run it with Python 3 and `matplotlib` to
regenerate its image in `img/`.

## Licence

The slides, their source and the figures are licensed under the Creative
Commons Attribution 4.0 International licence (CC BY 4.0,
<https://creativecommons.org/licenses/by/4.0/>): you may share and adapt them
for any purpose, provided you credit the author (Min Yang) and indicate any
changes. The bundled fonts keep their own free licences: Lato under the SIL
Open Font License 1.1, DejaVu Sans Mono under the DejaVu Fonts License.

The package itself is a separate work licensed under GPL-3.

## Citing

Yang, M., Biedermann, S. and Tang, E. (2013). On optimal designs for nonlinear
models: a general and efficient algorithm. *Journal of the American
Statistical Association* 108(504), 1411-1420.
<https://doi.org/10.1080/01621459.2013.806268>
