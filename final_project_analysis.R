library(tidyverse)
library(ncdf4)
library(lmtest)


# Read in NOAA tide gauge data (downloaded from their Water Levels Website) ---------------------------------------------

# Upload the NOAA tide guage data that was downloaded from https://tidesandcurrents.noaa.gov/stations.html?type=Water+Levels
msl_all <- read.csv("/Users/jordanasevigny/git/time-series-analysis-final-project/msl_noaa.csv")

# Switch Mean Sea Level from m from cm and subtract longterm mean to ensure we are comparing anomalies between all locations / data sources
msl_sd <- msl_all %>% filter(station=="9410170") %>% mutate(MSL=MSL*100)  %>% mutate(MSL = MSL-mean(MSL)) # m to cm and cacl anomaly
msl_mb <- msl_all %>% filter(station=="9413450") %>% mutate(MSL=MSL*100)  %>% mutate(MSL = MSL-mean(MSL)) # m to cm and cacl anomaly


# Plot  the msls
ggplot(msl_sd, aes(x=time, y=MSL)) +
  geom_line(linewidth = 0.6)
ggplot(msl_mb, aes(x=time, y=MSL)) +
  geom_line(linewidth = 0.6)


# Make into timeseries
# San Diego
sd_ts = ts(data = msl_sd$MSL, frequency = 12, start = c(1994, 1))
# Monterey
mb_ts = ts(data = msl_mb$MSL, frequency = 12, start = c(1994, 1))



# Read in model (UCSC ERA5) data -------------------------------------------------------

start_year = 1994
end_year = 2019

# Monterey Bay
mb_sea_level <- nc_open("/Users/jordanasevigny/git/time-series-analysis-final-projectzeta_timeseries_wcofs2_1994_2020_253_434.nc")

mb_zeta  <- ncvar_get(mb_sea_level, "zeta")
mb_ot  <- ncvar_get(mb_sea_level, "ocean_time")

# Make seconds into date then switch sea level to cm from m
origin <- as.POSIXct("1940-01-01 00:00:00", tz = "UTC")
dates <- as.Date(origin + mb_ot)
mb_df <- data.frame(
  date  = dates,
  year  = as.integer(format(dates, "%Y")),
  month = as.integer(format(dates, "%m")),
  day   = as.integer(format(dates, "%d")),
  zeta = mb_zeta*100 # m to cm
)

# Make monthly resolution then subtract longterm mean to ensure we are comparing anomalies between all locations / data sources
mb_mm <- mb_df %>%
  filter(year >= start_year & year <= end_year) %>%
  group_by(year, month) %>%
  summarise(mean(zeta)) %>%
  ungroup() %>%
  rename(monthly_mean='mean(zeta)') %>%
  mutate(monthly_mean = monthly_mean-mean(monthly_mean)) # Make anomoly

# Make timeseries
mb_ts_wcofs = ts(data = mb_mm$monthly_mean, frequency = 12, start = c(start_year, 1))

# San Diego
sd_sea_level <- nc_open("/Users/jordanasevigny/git/time-series-analysis-final-project/zeta_timeseries_wcofs2_1994_2020_314_290.nc")

sd_zeta  <- ncvar_get(sd_sea_level, "zeta")
sd_ot  <- ncvar_get(sd_sea_level, "ocean_time")

# Make seconds into date then switch sea level to cm from m
origin <- as.POSIXct("1940-01-01 00:00:00", tz = "UTC")
dates <- as.Date(origin + sd_ot)
sd_df <- data.frame(
  date  = dates,
  year  = as.integer(format(dates, "%Y")),
  month = as.integer(format(dates, "%m")),
  day   = as.integer(format(dates, "%d")),
  zeta = sd_zeta*100 # m to cm
)

# Make monthly resolution then subtract longterm mean to ensure we are comparing anomalies between all locations / data sources
sd_mm <- sd_df %>%
  filter(year >= start_year & year <= end_year) %>%
  group_by(year, month) %>%
  summarise(mean(zeta)) %>%
  ungroup() %>%
  rename(monthly_mean ='mean(zeta)') %>%
  mutate(monthly_mean = monthly_mean-mean(monthly_mean)) # Make anomoly

# Make timeseries
sd_ts_wcofs = ts(data = sd_mm$monthly_mean, frequency = 12, start = c(start_year, 1))



# Plotting timeseries -----------------------------------------------------

# MB vs SD, NOAA vs UCSC ERA 5
par(mfrow = c(2, 1))

plot.ts(mb_ts, main = "Timeseries of sea level at Monterey Bay by NOAA and UCSC ERA5", ylab="Sea level (cm)", xlab="Year", cex.axis = 1.5, cex.lab = 1.6, cex.main = 2)
lines(mb_ts_wcofs, col="red")
plot.ts(sd_ts, main = "Timeseries of sea level at San Diego by NOAA and UCSC ERA5", ylab="Sea level (cm)", xlab="Year", cex.axis = 1.5, cex.lab = 1.6, cex.main = 2)
lines(sd_ts_wcofs, col="red")

# MB - SD, NOAA vs UCSC ERA 5
par(mfrow = c(2, 1))

plot.ts(mb_ts-sd_ts, main = "Timeseries of sea level difference between Monterey Bay and San Diego (NOAA)", ylab="Sea level difference (cm)", xlab="Year")
plot.ts(mb_ts_wcofs-sd_ts_wcofs, main = "Timeseries of sea level difference between Monterey Bay and San Diego (UCSC ERA5)", ylab="Sea level difference (cm)", xlab="Year")

dev.off()

# Make difference dataframes
ts_diff <- mb_ts-sd_ts
ts_wcofs_diff <- mb_ts_wcofs-sd_ts_wcofs


# Decomposition -----------------------

# Decompose
y_stl = stl(ts_diff,s.window=7,t.window=25) # also tried 25
plot(y_stl, main = "NOAA Sea Level Difference (MB-SD) Decomposition")
hist(y_stl$time.series[, "remainder"])

y_stl_wcofs = stl(ts_wcofs_diff,s.window=7,t.window=25) # also tried 25
plot(y_stl_wcofs, main = "UCSC ERA5 Sea Level Difference (MB-SD) Decomposition")
hist(y_stl_wcofs$time.series[, "remainder"])

# Overlay decompositions
par(mfrow=c(4,1), mar=c(4,5,3,2))

years <- seq(start_year, end_year, length.out=312)

# raw
plot(years, ts_diff, type="l", col="orange", lwd=2, ylab="Raw (cm)", xlab="", main="Decomposition of NOAA and UCSC ERA5 Sea Level Differences (MB-SD)", cex.axis = 2, cex.lab = 2.5, cex.main = 3)
lines(years, ts_wcofs_diff, col="purple", lwd=2)
legend("topright",
       legend=c("NOAA", "UCSC ERA5"),
       col=c("orange","purple"),
       lwd=2,
       bty="n")

# seasonal
plot(years, y_stl$time.series[, "seasonal"], type="l", col="orange", lwd=2, ylab="Seasonal (cm)", xlab="", cex.axis = 2, cex.lab = 2.5, cex.main = 3)
lines(years, y_stl_wcofs$time.series[, "seasonal"], col="purple", lwd=2)
legend("topright",
       legend=c("NOAA", "UCSC ERA5"),
       col=c("orange","purple"),
       lwd=2,
       bty="n")

# trend
plot(years, y_stl$time.series[, "trend"], type="l", col="orange", lwd=2, ylab="Trend (cm)", xlab="", cex.axis = 2, cex.lab = 2.5, cex.main = 3)
lines(years, y_stl_wcofs$time.series[, "trend"],
      col="purple", lwd=2)
legend("topright",
       legend=c("NOAA", "UCSC ERA5"),
       col=c("orange","purple"),
       lwd=2,
       bty="n")

# remainder
plot(years, y_stl$time.series[, "remainder"], type="l", col="orange", lwd=2, ylab="Residual (cm)", xlab="Year", cex.axis = 2, cex.lab = 2.5, cex.main = 3)
lines(years, y_stl_wcofs$time.series[, "remainder"],
      col="purple", lwd=2)

legend("topright",
       legend=c("NOAA", "UCSC ERA5"),
       col=c("orange","purple"),
       lwd=2,
       bty="n")

# Plot remainder histograms
par(mfrow = c(1,2), mar=c(4,5,4,2))
hist(y_stl$time.series[, "remainder"], main = "NOAA Sea Level Difference (MB-SD)\nResidual Histogram", xlab="Residual (cm)", col= "orange", cex.axis = 2, cex.lab = 2.5, cex.main = 2)
hist(y_stl_wcofs$time.series[, "remainder"], main = "UCSC ERA5 Sea Level Difference (MB-SD)\nResidual Histogram", xlab="Residual (cm)", col="purple", cex.axis = 2, cex.lab = 2.5, cex.main = 2)
dev.off()

# Decomposition ranges

# Raw
max(y_stl$time.series[, "seasonal"]) - min(y_stl$time.series[, "seasonal"])
max(y_stl_wcofs$time.series[, "seasonal"]) - min(y_stl_wcofs$time.series[, "seasonal"])
quantile(y_stl$time.series[, "seasonal"])
quantile(y_stl_wcofs$time.series[, "seasonal"])

# Trend
max(y_stl$time.series[, "trend"]) - min(y_stl$time.series[, "trend"])
max(y_stl_wcofs$time.series[, "trend"]) - min(y_stl_wcofs$time.series[, "trend"])
quantile(y_stl$time.series[, "trend"])
quantile(y_stl_wcofs$time.series[, "trend"])

# Seasonal
max(ts_diff) - min(ts_diff)
max(ts_wcofs_diff) - min(ts_wcofs_diff)
quantile(ts_diff)
quantile(ts_wcofs_diff)

# Residual
max(y_stl$time.series[, "remainder"]) - min(y_stl$time.series[, "remainder"])
max(y_stl_wcofs$time.series[, "remainder"]) - min(y_stl_wcofs$time.series[, "remainder"])
quantile(y_stl$time.series[, "remainder"])
quantile(y_stl_wcofs$time.series[, "remainder"])

# Quartile table
tbl <- data.frame(
  source = c(rep("NOAA Sea Level Difference", 4), rep("UCSC ERA5 Sea Level Difference", 4)),
  element = c("raw", "trend", "seasonal", "residual", "raw", "trend", "seasonal", "residual"),
  min = c(quantile(ts_diff)[1], quantile(y_stl$time.series[, "trend"])[1], quantile(y_stl$time.series[, "seasonal"])[1], quantile(y_stl$time.series[, "remainder"])[1], quantile(ts_wcofs_diff)[1], quantile(y_stl_wcofs$time.series[, "trend"])[1], quantile(y_stl_wcofs$time.series[, "seasonal"])[1], quantile(y_stl_wcofs$time.series[, "remainder"])[1]),
  q25 = c(quantile(ts_diff)[2], quantile(y_stl$time.series[, "trend"])[2], quantile(y_stl$time.series[, "seasonal"])[2], quantile(y_stl$time.series[, "remainder"])[2], quantile(ts_wcofs_diff)[2], quantile(y_stl_wcofs$time.series[, "trend"])[2], quantile(y_stl_wcofs$time.series[, "seasonal"])[2], quantile(y_stl_wcofs$time.series[, "remainder"])[2]),
  median = c(quantile(ts_diff)[3], quantile(y_stl$time.series[, "trend"])[3], quantile(y_stl$time.series[, "seasonal"])[3], quantile(y_stl$time.series[, "remainder"])[3], quantile(ts_wcofs_diff)[3], quantile(y_stl_wcofs$time.series[, "trend"])[3], quantile(y_stl_wcofs$time.series[, "seasonal"])[3], quantile(y_stl_wcofs$time.series[, "remainder"])[3]),
  q75 = c(quantile(ts_diff)[4], quantile(y_stl$time.series[, "trend"])[4], quantile(y_stl$time.series[, "seasonal"])[4], quantile(y_stl$time.series[, "remainder"])[4], quantile(ts_wcofs_diff)[4], quantile(y_stl_wcofs$time.series[, "trend"])[4], quantile(y_stl_wcofs$time.series[, "seasonal"])[4], quantile(y_stl_wcofs$time.series[, "remainder"])[4]),
  max = c(quantile(ts_diff)[5], quantile(y_stl$time.series[, "trend"])[5], quantile(y_stl$time.series[, "seasonal"])[5], quantile(y_stl$time.series[, "remainder"])[5], quantile(ts_wcofs_diff)[5], quantile(y_stl_wcofs$time.series[, "trend"])[5], quantile(y_stl_wcofs$time.series[, "seasonal"])[5], quantile(y_stl_wcofs$time.series[, "remainder"])[5])
)

print(tbl)
write.csv(tbl, "/Users/jordanasevigny/git/time-series-analysis-final-project/Final project/quantile_table.csv")


# Plot UCSC ERA 5 vs NOAA for trend, seasonality, residual
# Raw
plot(ts_wcofs_diff, ts_diff)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
# Trend
plot(seq(1,312, by=1), y_stl$time.series[, "trend"])
plot(seq(1,312, by=1), y_stl_wcofs$time.series[, "trend"])
plot(y_stl_wcofs$time.series[, "trend"], y_stl$time.series[, "trend"])
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
# Seasonality
plot(y_stl_wcofs$time.series[, "seasonal"], y_stl$time.series[, "seasonal"])
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)


# Remainder / Interannual linear model and normal, independence, and variance tests
model = lm(y_stl$time.series[, "remainder"] ~ y_stl_wcofs$time.series[, "remainder"])
summary(model)
# 1.19 observation difference / model difference
# adj R2 = 0.75
# p-value: < 2.2e-16

# Remainder
plot(y_stl_wcofs$time.series[, "remainder"], y_stl$time.series[, "remainder"], xlab="NOAA Remainder", ylab="UCSC ERA5 Remainder", main="UCSC ERA5 vs NOAA MB-SD Sea Level Difference Interannual Variability (Residual)")
abline(a = 0, b = 1, col = "gray", lty = 2, lwd = 2)
abline(model, col="red", lwd=2)   
legend(
  "topleft",
  legend = c("1:1 line", "Linear regression"),
  col = c("gray", "red"),
  lty = c(2, 1),
  lwd = c(2, 2),
  bty = "n"
)

# Test assumptions
dwtest(model)
shapiro.test(model$residuals) 
bartlett.test(model$residuals, c(rep(0, length(model$residuals)/2), rep(1, length(model$residuals)/2))) 


# Cochrane-Orcutt transformation
y = y_stl$time.series[, "remainder"]
n = length(y)
rhohat = sum(model$residuals[-1]*model$residuals[-n])/sum(model$residuals^2)

# transform
ys = y[-1]-rhohat*y[-n]

# Fit a new model on the transformed data
model.t = lm(ys~y_stl_wcofs$time.series[, "remainder"][-1])
summary(model.t) 
# 1.18 observation difference / model difference
# p-value: < 2.2e-16

# Test assumptions using transformation
dwtest(model.t) # DW = 2.0101, p-value = 0.5326
shapiro.test(model.t$residuals) # W = 0.99707, p-value = 0.8485
bartlett.test(model.t$residuals, c(rep(0, 156), rep(1, 155))) # Bartlett's K-squared = 2.3915, df = 1, p-value = 0.122


# Remainder regression plot
par(cex.axis = 1.8, cex.lab = 2, mar=c(4,5,4,2))

plot(y_stl_wcofs$time.series[, "remainder"], y_stl$time.series[, "remainder"], xlab = "NOAA Residual (cm)", ylab = "UCSC ERA5 Residual (cm)", main = "UCSC ERA5 vs NOAA MB-SD Sea Level Difference\nInterannual Variability (Remainder)", cex.main = 2)

abline(a = 0, b = 1, col = "gray", lty = 2, lwd = 2)
abline(model.t, col="red", lwd=2)

legend(
  "topleft",
  legend = c("1:1 line", "Linear regression"),
  col = c("gray", "red"),
  lty = c(2, 1),
  lwd = c(2, 2),
  bty = "n"
)

# Plot residual timeseries 
par(mfrow = c(2, 1))
plot(time(mb_ts), y_stl$time.series[, "remainder"], type="l", ylab="Residual (cm)", xlab="Year", main="NOAA Sea Level Difference (MB-SD) Interannual Variability (remainder after detrending and deseasoning)")
abline(a=0, b=0, col="gray")
plot(time(mb_ts), y_stl_wcofs$time.series[, "remainder"], type="l", ylab="Residual (cm)", xlab="Year", main="UCSC ERA5 Sea Level Difference (MB-SD) Interannual Variability (remainder after detrending and deseasoning)")
abline(a=0, b=0, col="gray")


# Citations
citation()
citation("tidyverse")
citation("ncdf4")
citation("lmtest")
