# Introduction to population genetics — Activity 1

R code and data for the MBiol practical *How can spatial pattern of genetic diversity
inform conservation?*, which uses microsatellite data from Siamese rosewood
(*Dalbergia cochinchinensis*) across Indochina.

Rosewood is the most trafficked wild product in the world, and restoration planting is
part of the response. But a planting programme that ignores standing genetic diversity
can build a genetic bottleneck into the new forest. The practical asks where the seed
should come from, and makes students answer it from the data.

It is a two-hour self-led session in seven stations. Students load and check the data,
test Hardy–Weinberg and linkage, measure diversity corrected for sampling effort, then
describe population structure with pairwise *F*<sub>ST</sub>, a Mantel test and DAPC.
It ends with a conservation recommendation they have to defend.

This repository holds the code and the data. The student-facing handout is kept
separately.

## Contents

| File | What it is |
| --- | --- |
| `rosewood_activity1.Rmd` | The practical itself: every code block from the handout, Stations 0 to 6 |
| `rosewood_activity1_aux.Rmd` | Tutor's file: checks every result quoted in the handout, and rebuilds the teaching figures |
| `rosewood.csv` | Genotypes, GenAlEx format |
| `rosewood_coords.csv` | Longitude and latitude of each population |
| `rosewood_water.csv` | Rivers and lakes, for the map in Station 5 |

Both documents render to GitHub-flavoured markdown, so you can read a full worked run
without installing anything:

- **[rosewood_activity1.md](rosewood_activity1.md)** — the practical, with every output
- **[rosewood_activity1_aux.md](rosewood_activity1_aux.md)** — the checks and the figures

Those two files are generated. Edit the `.Rmd` and re-knit; do not edit the `.md` by
hand. The figures they refer to live in `outputs/` and `figures/`, and are committed so
the rendered pages display on GitHub.

## Requirements

R, and nine packages from CRAN:

```r
install.packages(c("poppr", "pegas", "hierfstat", "adegenet", "ade4",
                   "lattice", "mapplots", "maps", "RColorBrewer"))
```

`poppr` pulls in `adegenet`, but listing it explicitly keeps the setup script honest
about what the practical uses. Knitting the documents also needs `knitr` and
`rmarkdown`; RStudio supplies pandoc.

Built and checked with R 4.5.1 and poppr 2.9.8, pegas 1.3, hierfstat 0.5.11,
adegenet 2.1.11, ade4 1.7.23, lattice 0.22.7, mapplots 1.5.3, maps 3.4.3,
RColorBrewer 1.1.3, knitr 1.50, rmarkdown 2.29.

## Running it

Clone the repository, open either `.Rmd` in RStudio and click Knit. Both documents read
the CSVs by relative path, so nothing needs editing as long as the files stay in one
folder. Running them line by line works just as well, and is closer to what students
actually do.

Keep the clone somewhere with a short path. Knitting writes figures into
`rosewood_activity1_files/figure-gfm/`, which is long enough that a deeply nested
folder can trip the 260-character limit on Windows.

The practical takes about two minutes to knit; the auxiliary file takes longer, because
it repeats the whole analysis before drawing anything.

### One thing to know about Station 5.3

`find.clusters()` draws the BIC graph and then waits for you to type the number of
clusters. That is the point of the task: students read the graph and choose. When the
document is knitted there is nobody at the keyboard, so the chunk students run is shown
but not executed, and an identical chunk below it answers `5` on its own. The two give
the same clusters.

## The auxiliary file

`rosewood_activity1_aux.Rmd` is for whoever is running the session, not for students. It
does two things the practical deliberately leaves out:

- **Checks.** It re-runs the analysis and compares every result with the value quoted in
  the handout's *Expect to see* boxes, printing `PASS` or `FAIL` for each. There are 123
  of them. A failure almost always means a different package version.
- **Teaching figures.** It rebuilds the seven figures used in the handout and the
  pre-practical slides: the Wahlund effect, linkage disequilibrium, allele accumulation,
  *F*<sub>IS</sub> confidence intervals, isolation by distance, the BIC curve and a
  Mantel null distribution.

Running it after a package update is the quickest way to find out whether the numbers in
the handout still hold.

## The data

**`rosewood.csv`** — 523 trees, 9 microsatellite loci, 26 populations, from Hartvig et
al. In GenAlEx layout, so the first three rows are metadata: row 1 gives the number of
loci, individuals and populations followed by the sample sizes, row 2 the population
names, row 3 the column headings. Trees start at row 4, two columns per locus.

**`rosewood_coords.csv`** — one row per population: name, longitude, latitude. The
population names match `rosewood.csv` exactly, which matters, because Stations 4 and 5
join tables by name rather than by row order and the practical makes a point of it.

**`rosewood_water.csv`** — the Mekong and the other main rivers, Tonle Sap and three
large reservoirs, clipped from the Natural Earth 1:10m physical vectors (public domain)
to 100–111°E, 9–20°N. It is stored as plain longitude and latitude with a blank row
between features, so base R can draw it with `lines()` and `polygon()`. That keeps a
spatial dependency such as `sf` out of the practical. Columns are `type` (`river` or
`lake`), `name`, `lon` and `lat`.

## Reproducibility

Every step that samples at random is preceded by `set.seed(2026)`, so the numbers in the
handout are the numbers students get. The auxiliary file is what enforces that: if a
result drifts, a check fails and says by how much.

## Source

Hartvig I, So T, Changtragoon S, et al. (2018) Population genetic structure of the
endemic rosewoods *Dalbergia cochinchinensis* and *D. oliveri* at a regional scale
reflects the Indochinese landscape and life-history traits. *Ecology and Evolution*
8(1): 530–545. https://doi.org/10.1002/ece3.3626

With thanks to Prof Ida Hartvig (Copenhagen) and collaborators on the rosewood project,
whose data make this practical possible.
