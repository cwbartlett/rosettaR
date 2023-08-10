# rosettaR

`rosettaR` provides a method to combine datasets that measure the same latent traits when there is only partial overlap of measurements across the constituent datasets.  In this version it is possible to have measures that are local to only one dataset, though performance will be better if local variables are kept to a minimum.
## Installation

```{r}
devtools::install_github("cwbartlett/rosettaR")
```

## Examples

```{r}
library(rosetta)

# simulate data
d <- sim(seed = 100)

# check feature names
lapply(d$missing, names)


# run rosetta
d_rosetta <- rosetta(
  d = d$missing,
  factor_structure = list(
    a = c("a_1", "a_2", "a_3"),
    b = c("b_1", "b_2", "b_3"),
    c = c("c_1", "c_2", "c_3")
  )
)
```
