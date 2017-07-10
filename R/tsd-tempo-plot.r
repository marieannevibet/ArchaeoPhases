library(ggplot2)
library(ggthemes)
tsd.tempo.plot <- function (data, position, plot.result = NULL, level = 0.95,
                            count = TRUE, Gauss = FALSE, title = "Tempo plot",
                            subtitle = NULL, caption = NULL,
                            legend.title = "Legend",
                            legend.labels = c("Bayes estimate",
                                              "Credible interval, low",
                                              "Credible interval, high",
                                              "Gaussian approx., high",
                                              "Gaussian approx., low"),
                            x.label = "Calendar year",
                            y.label = "Cumulative events",
                            line.types = c("solid", "12", "11", "28", "28"),
                            plot.wi = 7, plot.ht = 7, plot.units = "in",
                            x.min = NULL, x.max = NULL, colors = TRUE,
                            out.file = NULL)
{
  if (is.null(plot.result))
    {
      L = length(position)
      groupOfDates = matrix(ncol = L, nrow = nrow(data))
      for (i in 1:L) {
          groupOfDates[, i] = data[, position[i]]
      }
      min = min(apply(groupOfDates, 2, min))
      max = max(apply(groupOfDates, 2, max))
      x = 10^c(0:10)
      if (min != 0) {
          c = 0
          for (i in 1:length(x)) {
              if (abs(min/x[i]) > 1) {
                  c = c + 1
              }
          }
          if (c > 3) {
              min = floor(min/x[c - 1]) * x[c - 1]
          }
          else {
              min = floor(min/x[c]) * x[c]
          }
      }
      if (max != 0) {
          d = 0
          for (i in 1:length(x)) {
              if (abs(max/x[i]) > 1) {
                  d = d + 1
              }
          }
          if (d > 3) {
              max = ceiling(max/x[d - 1]) * x[d - 1]
          }
          else {
              max = ceiling(max/x[d]) * x[d]
          }
      }
      t = seq(min, max, length.out = 50 * ncol(groupOfDates))
      f = function(x) {
          g = ecdf(x)
          y = g(t)
          if (count)
              y = y * ncol(groupOfDates)
          y
      }
      F = t(apply(groupOfDates, 1, f))
      moy = apply(F, 2, mean)
      ec = apply(F, 2, sd)
      qu = cbind(apply(F, 2, quantile, probs = (1 - level)/2, type = 8),
          apply(F, 2, quantile, probs = 1 - ((1 - level)/2), type = 8))
      quG = cbind(moy + qnorm(1 - (1 - level)/2) * ec, moy - qnorm(1 -
          (1 - level)/2) * ec)
      result = list(t = t, moy = moy, qu = qu, quG = quG)
    }
  else
  {
    result = plot.result
  }
      if (Gauss) {
          result.mat <- cbind(result$moy, result$qu, result$quG)
          colnames(result.mat) <- legend.labels
      }
      else {
          result.mat <- cbind(result$moy, result$qu)
          colnames(result.mat) <- legend.labels[1:3]
      }
      plot.result <- as.data.frame.table(result.mat)
      colnames(plot.result) <- c("Var1", "Legend", "Count")
      plot.result$Year <- result$t
      if (colors) {
          h <- ggplot(plot.result, aes(x = plot.result$Year, y = plot.result$Count,
              colour = plot.result$Legend))
          h <- h + guides(colour = guide_legend(title = legend.title))
      }
      else {
          h <- ggplot(plot.result, aes(x = plot.result$Year, y = plot.result$Count,
              linetype = plot.result$Legend))
          h <- h + scale_linetype_manual(values = line.types,
                                         guide = guide_legend(title = legend.title))
      }
      h <- h + geom_line()
      h <- h + scale_y_continuous(breaks = pretty(x = plot.result$Count))
      h <- h + labs(x = x.label,
                    y = y.label,
                    title = title,
                    subtitle = subtitle,
                    caption = caption)
      if (!is.null(x.min) & !is.null(x.max)) {
        h <- h + xlim(x.min, x.max)
      }
      if (!is.null(out.file)) {
          ggsave(filename = out.file, plot = h, height = plot.ht,
              width = plot.wi, units = plot.units)
      }
      dev.new(height = plot.ht, width = plot.wi)
      print(h)
      result
  }
## theme_set(theme_bw())
theme_set(theme_solarized_2(light = FALSE))
scale_colour_discrete <- scale_colour_solarized
res <- tsd.tempo.plot(habitation.events, habitation.index, plot.result = res,
               title = "Habitation Construction", colors = TRUE, plot.ht = 3.5,
               plot.wi = 7.8, caption = "version 1.0", x.min = 1000, x.max = 2000,
               out.file = "~/Public/projects/949-current-anth/figure/habitation-test.pdf")

res <- tsd.tempo.plot(data, c(2:12), plot.result = NULL,
                      title = "TempoPlot Cuba", colors = TRUE)#, plot.ht = 3.5,
                      plot.wi = 7.8, caption = "version 1.0")#, x.min = 1000, x.max = 2000)
#,
                     # out.file = "~/Public/projects/949-current-anth/figure/habitation-test.pdf")



