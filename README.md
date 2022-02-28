# nexonsAnalysis

R functions for analysis and graphical representation of output from
nexons

``` r
devtools::install_github("laurabiggins/nexonsAnalysis")
library(nexonsAnalysis)
```

# Example data and plot

## Using nexons gtf output

Nexons can produce output in 2 different formats. One we will refer to
as the “default nexons output”, the other is in gtf format. This package
was originally designed to work with the gtf format, but it can also use
the default nexons output. The example below using the `plot_wrapper()`
function is a very simple example that does not allow many arguments to
be specified and will only work with the nexons gtf file output.

``` r
file <- system.file("extdata", "nexons_sirv5_f15.gtf", package = "nexonsAnalysis")
plot_wrapper(file, min_count = 5, quant_plot = TRUE)
```

    FALSE [1] "No ordering of y axis specified"

![](man/figures/README-unnamed-chunk-2-1.png)

The above example can be reproduced using the individual functions that
`plot_wrapper()` is a wrapper around. This allows more options to be
specified.

``` r
file <- system.file("extdata", "nexons_sirv5_f15.gtf", package = "nexonsAnalysis")
nexons_output <- readr::read_delim(file)
parsed_splices <- parse_nexons_gtf(nexons_output, min_count = 3)
draw_splice_picture(parsed_splices, quant = TRUE, order_splices = "score", gene="SIRV5")
```

![](man/figures/README-unnamed-chunk-3-1.png)

## Using default nexons output

### For one sample at a time

``` r
file <- system.file("extdata", "sirv5.txt", package = "nexonsAnalysis")
nexons_output <- readr::read_delim(file)
# To plot the splice pictures for 1 sample
parsed_splices <- parse_default_nexons(nexons_output, score_column = "seqs_sirv5_minimap.sam")
draw_splice_picture(parsed_splices, quant = TRUE, order_splices = "name")
```

![](man/figures/README-unnamed-chunk-4-1.png)

### For multiple samples

Use purrr::map to iterate over multiple datasets.

``` r
file <- system.file("extdata", "sirv5.txt", package = "nexonsAnalysis")
nexons_output <- readr::read_delim(file)
# get names of datasets - there are 2 in this file
count_columns <- tail(colnames(nexons_output), n=2)
parsed_splices <- purrr::map(count_columns, parse_default_nexons, nexons_output=nexons_output)
purrr::map(parsed_splices, draw_splice_picture, quant=TRUE, order_splices = "score")
```

![](man/figures/README-unnamed-chunk-5-1.png)![](man/figures/README-unnamed-chunk-5-2.png)

    FALSE [[1]]
    FALSE TableGrob (1 x 2) "arrange": 2 grobs
    FALSE   z     cells    name           grob
    FALSE 1 1 (1-1,1-1) arrange gtable[layout]
    FALSE 2 2 (1-1,2-2) arrange gtable[layout]
    FALSE 
    FALSE [[2]]
    FALSE TableGrob (1 x 2) "arrange": 2 grobs
    FALSE   z     cells    name           grob
    FALSE 1 1 (1-1,1-1) arrange gtable[layout]
    FALSE 2 2 (1-1,2-2) arrange gtable[layout]
