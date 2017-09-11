## ---- echo = FALSE, message = FALSE-----------------------------------------------------------------------------------
knitr::opts_chunk$set(comment = "")
options(width = 120, max.print = 100)
library(ArchaeoPhases)

## ----fig.align='center',fig.width=6,fig.height=5----------------------------------------------------------------------
data("KADatesChronoModel")
TempoPlot(KADatesChronoModel, c(2:17), level = 0.95)

## ----fig.align='center',fig.width=6,fig.height=5----------------------------------------------------------------------
data("KADatesChronoModel")
TempoActivityPlot(KADatesChronoModel, c(2:17), level = 0.95)

## ----fig.align='center',fig.width=6,fig.height=5----------------------------------------------------------------------
data("KADatesChronoModel")
OccurrencePlot(KADatesChronoModel, c(2:17), level = 0.95, print.data.result = FALSE)

