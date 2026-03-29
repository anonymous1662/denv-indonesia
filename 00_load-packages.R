# code/00_load-packages.R

library(here)
library(dplyr)
library(janitor)
library(tidyr)
library(readxl)
#library(ggplot2)
#library(showtext)
library(stringr)
library(tibble)
library(yaml)
library(qs)
library(lubridate)
library(abind)
library(data.table)
#library(igraph)
#library(future.apply)
#library(float)
library(odin2)
library(dust2)
#library(decor)
#library(pkgbuild)
#library(pkgload)

# --- Plot theme ---
#theme_set(theme_bw() +
#            theme(
#              plot.title        = element_text(size = 8, family = "sans"),
#              axis.title        = element_text(size = 8, family = "sans"),
#              axis.text         = element_text(size = 7, family = "sans"),
#              legend.title      = element_text(size = 8, family = "sans"),
#              plot.subtitle     = element_text(size = 8, family = "sans"),
#              legend.text       = element_text(size = 7, family = "sans"),
#              legend.key.height = unit(0.5, "cm"),
#              legend.position   = "bottom",
#              strip.background  = element_rect(fill = "#082544"),
#              strip.text        = element_text(color = "white", size = 6, family = "sans"),
#              panel.border      = element_rect(color = "#082544", fill = NA, linewidth = 0.5),
#              panel.grid.major  = element_line(color = "grey90", linewidth = 0.3),
#              panel.grid.minor  = element_line(color = "grey95", linewidth = 0.2),
#              axis.line         = element_blank(),
#              axis.ticks        = element_line(color = "#082544")
#            )
#)



