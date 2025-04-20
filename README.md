# Gene Isoform Visualizer App

This is a self-study project to practice programming in R, specifically using the **Shiny** package.

## Motivation

Whenever I needed to download data from Ensembl, I often felt confused by the many filters and options in BioMart. Creating a tool that can retrieve data quickly and intuitively seemed like a good idea—especially for researchers or research groups that frequently work with a specific type of data.

## What the App Does

- You type in a gene name (e.g., `TP53`)
- The app queries **Ensembl BioMart** for the transcript and exon data using the `biomaRt` R package
- The alternative splicing isoforms are visualized using the `Gviz` package
- A downloadable table shows exon-level information, which can be saved as a `.csv` file

## Why it is not online

Unfortunately, I couldn't successfully deploy the app to [shinyapps.io](https://www.shinyapps.io). This is most likely due to the memory limitations of the free tier, which may be a problem when querying large datasets live from Ensembl. I also attempted to make the app use a pre-saved, full Ensembl dataset, but this approach didn't work either.

**For this reason, the app needs to be run locally.**

## How to Run Locally

Before usage, the required R packages should be installed:

```r
install.packages(c("shiny", "biomaRt", "Gviz", "DT", "dplyr", "shinycssloaders"))
