# T-cell receptor analysis of synucleinopathy

This repository holds scripts for generating figures for the manuscript:
Intestinal Macrophages Modulate Brain Synucleinopathy in Parkinson’s disease.

## Installation

To run the scripts, first install the local package via running the following
command in the top level directory:

```shell
pip install .
```

## Structure

All plots used in the manuscript are located within `./publication/`. Additional
plots are located in `./additional/`.

## Data

All data used in these plots can be obtained from the SRA: PRJNA1321765. To
obtain aligned repertoires, use
![decombinator](https://github.com/innate2adaptive/Decombinator).

Example command:

```shell
decombinator pipeline -c <chain> -in <filename> -br R2 -sp mouse -bl 42 -ol M13
```

See the readme of
![decombinator](https://github.com/innate2adaptive/Decombinator) for further
details.
