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
specified. The `order_splices` argument in draw_splice_picture accepts
one of (one of ‘score’, ‘name’, NULL) default NULL. This refers to the
ordering of the splices on the y axis (on the left hand plot), ‘score’
will sort data by score (highest at the top), ‘name’ will sort data
alphabetically (with unknowns at the bottom). `draw_splice_picture` also
allows a particular gene to be specified with the `gene` argument.

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
p <- purrr::map(parsed_splices, draw_splice_picture, quant=TRUE, order_splices = "score")
```

![](man/figures/README-unnamed-chunk-5-1.png)![](man/figures/README-unnamed-chunk-5-2.png)

### Adding titles to the plots

``` r
file <- system.file("extdata", "sirv5.txt", package = "nexonsAnalysis")
nexons_output <- readr::read_delim(file)
# get names of datasets - there are 2 in this file
count_columns <- tail(colnames(nexons_output), n=2)
parsed_splices <- purrr::map(count_columns, parse_default_nexons, nexons_output=nexons_output) |>
  purrr::set_names(count_columns)
p <- purrr::imap(parsed_splices, ~ draw_splice_picture(.x, title_text=.y, quant=TRUE, order_splices = "score"))
```

![](man/figures/README-unnamed-chunk-6-1.png)![](man/figures/README-unnamed-chunk-6-2.png)

## Flagging up potentially truncated reads

The function `identifyPotentialTruncations` identifies whether
transcripts might be truncations of other transcripts. It adds an
additional column to the tibble named truncation_origin.

``` r
file <- system.file("extdata", "nexons_sirv5_f15_trunc.txt", package = "nexonsAnalysis")
nexons_output <- readr::read_delim(file)
parsed_splices <- parse_nexons_gtf(nexons_output, min_count = 3)

parsed_with_trunc <- identifyPotentialTruncations(parsed_splices, flexibility = 10)
knitr::kable(parsed_with_trunc)
```

| strand | score | Transcript_id | Gene_id | splice_pattern                                                                                                                                                               | variant | truncation_origin             |
|:--|--:|:----|:---|:----------------------------------------------|---:|:--------|
| \+     |    99 | SIRV501       | SIRV5   | 1149:1988-2033:2120-2315:3299-3404:3484-3643:5381-5450:5544-5626:6112-6169:6328-6452:6659-6722:6827-6957:7145-7307:7682-7762:7871-8016:8278-8381:8455-8585:10859             |       1 |                               |
| \+     |    96 | SIRV502       | SIRV5   | 1149:1988-2033:2120-2156:2271-2488:3299-3404:3484-3643:5381-5450:5544-5626:6112-6169:6328-6452:6659-6722:6827-6957:7145-7307:7682-7762:7871-8016:8278-8381:8455-8585:10859   |       2 |                               |
| \+     |   100 | SIRV504       | SIRV5   | 13606                                                                                                                                                                        |       3 |                               |
| \+     |    92 | SIRV505       | SIRV5   | 1149:1988-2033:2120-2156:2271-2315:3299-3404:3484-3643:5381-5450:5544-5626:6112-6169:6328-6452:6827-6957:7145-7307:7682-7762:7871-8381:8455-8585:10859                       |       4 |                               |
| \+     |   200 | SIRV506       | SIRV5   | 1149:1988                                                                                                                                                                    |       5 | 1, 2, 4, 7, 9, 12, 13, 15, 16 |
| \+     |    95 | SIRV507       | SIRV5   | 1149:1926-2033:2120-2156:2271-2315:3299-3404:3484                                                                                                                            |       6 |                               |
| \+     |    94 | SIRV508       | SIRV5   | 1149:1988-2033:2120-2156:2271-2315:3299-3404:3484-3643:5381-5450:5544-5626:6112-6169:6328-6452:6659-6722:6827-6957:7145-7307:7682-7762:7871-8381:8455-8585:10859             |       7 |                               |
| \+     |   100 | SIRV509       | SIRV5   | 8381:8455-8585:10859-10991:11312                                                                                                                                             |       8 |                               |
| \+     |    93 | SIRV510       | SIRV5   | 1149:1988-2033:2120-2156:2271-2315:3299-3404:3484-3643:5381-5450:5544-5626:6112-6169:6328-6452:6827-6957:7145-7307:7682-7762:7871-8016:8278-8381:8455-8585:10859-10991:11134 |       9 |                               |
| \+     |   100 | SIRV512       | SIRV5   | 2406                                                                                                                                                                         |      10 |                               |
| \+     |   100 | unknown       | SIRV5   | 8585:10859                                                                                                                                                                   |      11 | 1, 2, 4, 7, 8, 9, 15, 16      |
| \+     |    40 | unknown       | SIRV5   | 1154:1987-2040:2123-2156:2274-2315:3299-3409:3484-3643:5381-5450:5544                                                                                                        |      12 | 4, 7, 9                       |
| \+     |    23 | unknown       | SIRV5   | 1149:1988-2033:2120-2156:2271-2315:3299-3404:3484                                                                                                                            |      13 | 4, 7, 9, 12                   |
| \+     |     4 | unknown       | SIRV5   | 1149:1926-2047:2265-2315:3299-3404:3484                                                                                                                                      |      14 |                               |
| \+     |     3 | unknown       | SIRV5   | 1149:1986-2047:2265-2315:3299-3404:3484-3643:5381-5450:5544-5626:6112-6169:6328-6452:6827-6957:7145-7307:7682-7776:7881-8381:8455-8585:10859                                 |      15 |                               |
| \+     |     3 | unknown       | SIRV5   | 1149:1988-2047:2265-2315:3299-3404:3484-3643:5381-5450:5544-5626:6112-6169:6328-6452:6827-6957:7145-7307:7682-7762:7871-8016:8280-8381:8455-8585:10859-10991:11134           |      16 |                               |
