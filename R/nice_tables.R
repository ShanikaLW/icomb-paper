# A function to create summary table
# Arguments
# lst: a named list of 2D matrices where the named rows give level wise performances across different 
# forecast horizons (named columns) for different methods. 
# caption: A caption for the table
# digits: number of decimal places to display
# relative: compute prial relative to this 2D matrix
# fix_bold_minus: fix the size of numbers after making them bold and size of the minus sign

library(kableExtra)
library(tidyverse)

bold_minus <- function(x) 
  ifelse(as.numeric(x) == min(as.numeric(x)), paste0("$\\pmb{", x, "}$"), paste0("$", x, "$"))

default <- function(x) 
  cell_spec(x, bold = ifelse(as.numeric(x) == min(as.numeric(x)), TRUE, FALSE))

make_nice_tables <- function(lst, caption = "PRIAL from ARIMA models", digits = 2, 
                             relative = NULL, fix_bold_minus = TRUE) {
  summary_list <- lst
  if (!is.null(relative)) 
    summary_list <- map(summary_list, function(x, y) 100 * (x/y - 1), y = relative)
  
  summary_list <- map(summary_list, round, digits = digits)
  level_names <- rownames(summary_list[[1]])
  dims <- dim(summary_list[[1]])
  summary_levelwise <- lapply(1:dims[1], function(i, x) t(sapply(x, "[", i, TRUE)), x = summary_list)
  add_levels <- setNames(rep(length(summary_list), dims[1]), level_names)
  map(summary_levelwise, ~.x %>%
        as_tibble(rownames = "Methods") %>%
        mutate(across(-Methods, ~as.character(formatC(.x, digits = digits, format = 'f')))) %>% 
        mutate(across(-Methods, ifelse(fix_bold_minus, bold_minus, default)))) %>%
    bind_rows() %>% 
    kbl(booktabs = TRUE,
        caption = caption,
        linesep = "" ,
        align = c("l", rep("r", dims[2])),
        format.args = list(nsmall = 3),
        escape = FALSE) %>%
    kable_styling(latex_options = c("HOLD_position", "scale_down"), 
                  bootstrap_options = c("condensed"), position = "center") %>% 
    pack_rows(index = add_levels)
  
}

# Test cases
# x <- array(rnorm(24, 30), dim = c(4, 6))
# y <- array(rnorm(24, 50), dim = c(4, 6))
# z <- array(rnorm(24, 60), dim = c(4, 6))
# a <- array(rnorm(24, 80), dim = c(4, 6))
# colnames(x) <- colnames(y) <- colnames(z) <- colnames(a) <- c(1:5, "Average")
# rownames(x) <- rownames(y) <- rownames(z) <- rownames(a) <- c("Top", "Level 1", "Bottom", "Average")
# 
# make_nice_tables(BU = x, `WLS$_v$` = y, OLS = z)
# make_nice_tables(BU = x, `WLS$_v$` = y, OLS = z, relative = a)
