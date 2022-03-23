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

## 

Flagging up potentially truncated reads

``` r
devtools::load_all()
file <- system.file("extdata", "inst/extdata/nexons_sirv5_f15_trunc.txt", package = "nexonsAnalysis")
file <- "inst/extdata/nexons_sirv5_f15_trunc.txt" #!! temp - remove this
nexons_output <- readr::read_delim(file)
parsed_splices <- parse_nexons_gtf(nexons_output, min_count = 3)

(parsed_with_trunc <- identifyPotentialTruncations(parsed_splices))
```

    FALSE # A tibble: 16 x 7
    FALSE    strand score Transcript_id Gene_id splice_pattern    variant truncation_orig~
    FALSE    <chr>  <dbl> <chr>         <chr>   <chr>               <dbl> <chr>           
    FALSE  1 +         99 SIRV501       SIRV5   1149:1988-2033:2~       1 ""              
    FALSE  2 +         96 SIRV502       SIRV5   1149:1988-2033:2~       2 ""              
    FALSE  3 +        100 SIRV504       SIRV5   13606                   3 ""              
    FALSE  4 +         92 SIRV505       SIRV5   1149:1988-2033:2~       4 ""              
    FALSE  5 +        200 SIRV506       SIRV5   1149:1988               5 "1, 2, 4, 7, 9,~
    FALSE  6 +         95 SIRV507       SIRV5   1149:1926-2033:2~       6 ""              
    FALSE  7 +         94 SIRV508       SIRV5   1149:1988-2033:2~       7 ""              
    FALSE  8 +        100 SIRV509       SIRV5   8381:8455-8585:1~       8 ""              
    FALSE  9 +         93 SIRV510       SIRV5   1149:1988-2033:2~       9 ""              
    FALSE 10 +        100 SIRV512       SIRV5   2406                   10 ""              
    FALSE 11 +        100 unknown       SIRV5   8585:10859             11 "1, 2, 4, 7, 8,~
    FALSE 12 +         40 unknown       SIRV5   1154:1987-2040:2~      12 "4, 7, 9"       
    FALSE 13 +         23 unknown       SIRV5   1149:1988-2033:2~      13 "4, 7, 9, 12"   
    FALSE 14 +          4 unknown       SIRV5   1149:1926-2047:2~      14 ""              
    FALSE 15 +          3 unknown       SIRV5   1149:1986-2047:2~      15 ""              
    FALSE 16 +          3 unknown       SIRV5   1149:1988-2047:2~      16 ""
