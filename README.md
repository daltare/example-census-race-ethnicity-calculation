This repository contains source files for a document (published as a GitHub pages site [here](https://daltare.github.io/example-census-race-ethnicity-calculation)) that provides an example of how to use tools available from the [R programming language](https://www.R-project.org/) to estimate characteristics of any given *target* spatial area(s) (e.g., neighborhoods, project boundaries, water supplier service areas, etc.) based on data from a *source* dataset containing the characteristic data of interest (e.g., census data, CalEnvrioScreen scores, etc.), especially when the boundaries of the *source* and *target* areas overlap but don't necessarily align with each other. It also provides some brief background on the various types of data available from the U.S Census Bureau, and links to a few places to find more in-depth information.

### Publishing

This document is published as a site on GitHub pages, following the instructions [here](https://quarto.org/docs/publishing/github-pages.html). In particular:

-   Start [here](https://quarto.org/docs/publishing/github-pages.html#source-branch) to set up the `gh-pages` branch in Git / GitHub and format `.gitignore` to ignore rendered directories
-   Then follow instructions [here](https://quarto.org/docs/publishing/github-pages.html#github-action) to set up a GitHub Action for publication

More details are available [here](https://github.com/daltare/test-gh-actions#quarto).
