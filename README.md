
<!-- README.md is generated from README.Rmd. Please edit that file -->
carbon14
========

The goal of carbon14 is to provide a tidy interface to radiocarbon dating, to promote its inclusion in reproducible manuscripts.

Installation
------------

You can install carbon14 from github with:

``` r
# install.packages("devtools")
devtools::install_github("paleolimbot/carbon14")
```

Example
-------

Radiocarbon ages from [Long Lake, Nova Scotia, Canada](http://www.facetsjournal.com/doi/10.1139/facets-2017-0004):

    #> # A tibble: 5 x 4
    #>   sample_id       depth age_14C age_error
    #>   <chr>           <dbl>   <dbl>     <dbl>
    #> 1 LL082011-C2-39  94.0     9148      49.0
    #> 2 UOC-0844        70.5     8582      28.0
    #> 3 LL082011-C2-87  46.0     4396      55.0
    #> 4 UOC-0845        37.5      575      18.0
    #> 5 LL082011-C2-124  9.00     623      34.0

Calibrating dates using `calibrate()`:

``` r
result <- dates %>%
  calibrate(
    measured_age = age_14C, 
    measured_age_error = age_error,
    name = sample_id
  )
result
#> # A tibble: 5 x 9
#>   name   measured_age measured_age_err…    df curve   cal_age   curve_name
#>   <chr>         <dbl>             <dbl> <dbl> <S3: a> <S3: cdi> <chr>     
#> 1 LL082…         9148              49.0   Inf <age_c… 10318.48… intcal13  
#> 2 UOC-0…         8582              28.0   Inf <age_c… 9539.570… intcal13  
#> 3 LL082…         4396              55.0   Inf <age_c… 5005.665… intcal13  
#> 4 UOC-0…          575              18.0   Inf <age_c… 591.8846… intcal13  
#> 5 LL082…          623              34.0   Inf <age_c… 604.1533… intcal13  
#> # ... with 2 more variables: measured_age_type <chr>, cal_age_type <chr>
```

The result is a tibble, which also has a plot method to examine the results:

``` r
plot(result)
```

![](README-plot_long_lake-1.png)

And a summary of the ages can be generated using `summary()` on the `cal_age` column of the output:

``` r
summary(result$cal_age)
#> # A tibble: 5 x 8
#>   .name       weighted_mean quantile_5 quantile_25 quantile_50 quantile_75
#>   <chr>               <dbl>      <dbl>       <dbl>       <dbl>       <dbl>
#> 1 LL082011-C…         10318      10237       10263       10305       10362
#> 2 UOC-0844             9540       9528        9535        9539        9543
#> 3 LL082011-C…          5006       4870        4917        4975        5040
#> 4 UOC-0845              592        544         553         607         618
#> 5 LL082011-C…           604        555         575         599         635
#> # ... with 2 more variables: quantile_95 <dbl>, dist <chr>
```
