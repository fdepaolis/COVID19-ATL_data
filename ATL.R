# ## ## ## ## ## ## ## ##
# Project using The Tracking COVID Project data
# with API and JSON
#
#  Data by STATE
#
## ## ## ## ## ## ## ## ##


# The Atlantic  uses >>>>> "The COVID Tracking Project"
# https://github.com/COVID19Tracking 


# highlight:  should be one of “default”, “tango”, “pygments”, “kate”, “monochrome”, “espresso”, “zenburn”, “haddock”, “breezedark”, “textmate”

# New way of loading packages pacman::P-load()
# p_funs(i) gives you the functions available in package "i"
# pacman::p_path()  path to your R-library listing all packages

pacman::p_load(httr, jsonlite,dplyr, ggplot2, reshape2, ggplot2, forecast, fpp2)

# API
# territ <- c("Guam" , "Virgin Islands" , "Puerto Rico" , "Northern Mariana Islands", "American Samoa")

ATL_API_curr = GET("https://covidtracking.com/api/v1/states/current.json")  # usually a day or 2 behind

ATL_API_data_curr = fromJSON(rawToChar(ATL_API_curr$content))

names(ATL_API_data_curr)

ATL_API_hist = GET("https://covidtracking.com/api/v1/states/daily.json")  # usually a day or 2 behind

ATL_API_data_hist = fromJSON(rawToChar(ATL_API_hist$content))


# removes territories
territ <- c("AS","GU","MP","PR","VI")
ATL_API_data_curr_2 <- ATL_API_data_curr
# 
for (i in num){
ATL_API_data_curr_2 <- subset(ATL_API_data_curr_2,
        ATL_API_data_curr_2$state !=territ[i]
        )
}

## ## ## ## ## ## ## ## ##
# Projections based on linearization algorithms
## ## ## ## ## ## ## ## ##

CALIF_1 <- subset(ATL_API_data_hist, state=='CA', select=c(date,state,positive))
CALIF_2 <-CALIF_1 %>%
            arrange(date)   # sort ascending AND changes row number

# Transformations

CALIF_2$LN.y   <- log(CALIF_2$positive)
CALIF_2$days   <- seq.int(nrow(CALIF_2))-1

# Exponential model LN.y=f(x)
Exp.model_1 <- lm(CALIF_2$LN.y~CALIF_2$days)

# Compute the constant “a” for the exponential model
a.exponential<-exp(Exp.model_1$coeff[1]) 

#    (Intercept) 
#       812.0893  

# y = e^(ax); where a = 812.0893 and x is 'days'

plot(CALIF_2$positive)

# exponential fitting
CALIF_2[nrow(CALIF_2)+3,] <- NA

CALIF_2$Future.x <- seq.int(nrow(CALIF_2))-1

CALIF_2$y.exp<-a.exponential*exp(Exp.model_1$coeff[2]*CALIF_2$Future.x)

plot(CALIF_2$y.exp)


# Power
CALIF_2$y53     <- CALIF_2$positive - 53
CALIF_2$LN.y53  <- log(CALIF_2$y53)
CALIF_2$LN.x    <- log(CALIF_2$days)

Power.Model <- lm(CALIF_2$LN.y53[3:135] ~ CALIF_2$LN.x[3:135])
a.power<-exp(Power.Model$coeff[1]) 

#    (Intercept)  (coeff)
#      -0.15443  2.65079

# y = ax^b; x is 'days'; a = exp(intercept)=0.8569034 and b = 2.65079



CALIF_2$y.power <- a.power*(CALIF_2$Future.x^Power.Model$coeff[2])+CALIF_2$positive[1]


p <- ggplot() + 
  geom_line(data = CALIF_2, aes(days, positive/1000, color = "Actual"), size = 1.15) + 
  geom_line(data = CALIF_2, aes(days, y.power/1000, color = "Power fitting"), size = 0.9) +
  labs(title="California-Positive Cases",
          x ="Days since recording", 
          y = "Total Cases [1,000]") +
  scale_color_manual(name = "Cases", 
                     values = c("Actual" = "#999300", "Power fitting" = "#B4446B"))


p + theme(panel.background = element_rect(fill='gray90', colour='#999999'))

# COLOR SCHEME SELECTOR
# http://colorschemedesigner.com/csd-3.5/
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf 


# colors
# Mustard   "#E69F00"
# Deep gray "#999999"

library(colortools)
splitComp("steelblue")
tetradic("steelblue")
square("steelblue")
analogous("#00AF71")
complementary("#005496")

library("RColorBrewer")
display.brewer.all()


## ## ## ## ## ## ## ## ##
# Using package "fpp2"
## ## ## ## ## ## ## ## ##

library(fpp2)

# video https://www.youtube.com/watch?v=dBNy_A6Zpcc&t=717s

# Create several time-series for variable/place combination
# Start with 'historical data'
# California dataset

CAL_positiv <- CALIF_2[,3] # Extracts the column containing the number of positive cases

end_date <- paste(substr(as.character(CALIF_1[1, 1]), 1, 4),substr(as.character(CALIF_1[1, 1]), 5, 6),substr(as.character(CALIF_1[1, 1]), 7, 8), sep = "-") # Updates the latest date to use as end date in time series

inds <- seq(as.Date("2020-03-04"), as.Date(end_date), by = "day") # Creates a daily Date object - helps my work on dates

CAL_ts <- ts(CAL_positiv, 
           start = c(2020, as.numeric(format(inds[1], "%j"))),
           frequency = 365) # Creates a time series object


# detrend with first difference
CAL_1_diff <- diff(CAL_ts)
# log 1st diff
CAL_1_diff_log <- diff(log(CAL_ts))
# second difference
CAL_2_diff <- diff(diff(CAL_ts))
# log 2nd diff
CAL_2_diff_log <- diff(diff(log(CAL_ts)))
# third difference
CAL_3_diff <- (diff(CAL_2_diff))

# fitting
#naive forecast (snaive is for series with seasonality)
fit_CAL <- naive(CAL_1_diff)  
print(summary(fit_CAL))   # Residual sd: 1147.9273 
checkresiduals(fit_CAL)

## ## ## ## ##
# fit ETS model
## ## ## ## ##
fit_ets <- ets(CAL_ts)
print(summary(fit_ets))   # sigma:  1005.318 (equiv.Residual sd )
checkresiduals(fit_ets)  
# ETS(error, trend, seasonal) -> ETS(A,A,N) <=> "Holt’s linear method with additive errors" (from https://robjhyndman.com/eindhoven/2-1-StateSpace.pdf)

## ## ## ## ##
# ARIMA model
## ## ## ## ##
fit_arima <- auto.arima(CAL_ts,d=1, 
                        stepwise = FALSE, 
                        approximation = FALSE, 
                        trace = TRUE) 
# d=1 -> Best model: ARIMA(3,1,0) with drift  
# d=2 -> Best model: ARIMA(2,2,2)
# d=3 -> Best model: ARIMA(2,3,3)
print(summary(fit_arima))   # sigma^2 estimated as 1055447 or sigma = 1027.35
checkresiduals(fit_arima)       


fit_arima_log <- auto.arima(log(CAL_ts),d=1, 
                        stepwise = FALSE, 
                        approximation = FALSE, 
                        trace = TRUE) # Best model: ARIMA(5,1,0) with drift  
print(summary(fit_arima_log))   # sigma^2 estimated as 0.001426
checkresiduals(fit_arima_log)       

## ## ## ## ##
# Forecast
## ## ## ## ##
fcast <- forecast(fit_arima, h=7)
autoplot(fcast, include = 25)

print(summary(fcast))

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##

# Forecasts:
#             Point Forecast  Lo 80       Hi 80       Lo 95     Hi 95
# 2020.07-28    491448.9   486574.9    496322.9    483994.8  498903.0
# 2020.07-29    499017.5   492659.7    505375.3    489294.1  508740.9
# 2020.07-30    506434.1   498517.4    514350.8    494326.5  518541.6
# 2020.07-31    513805.8   504226.3    523385.3    499155.2  528456.3
# 2020.08-01    521148.0   509806.7    532489.2    503803.0  538493.0
# 2020.08-02    528417.6   515237.1    541598.0    508259.8  548575.4
# 2020.08-03    535628.5   520534.2    550722.8    512543.7  558713.3

## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##


#   ARIMA models for time series forecasting   
#   https://people.duke.edu/~rnau/411arim.htm
# A nonseasonal ARIMA model is classified as an "ARIMA(p,d,q)" model, where:
#         p is the number of autoregressive terms,
#         d is the number of nonseasonal differences needed for stationarity, and
#         q is the number of lagged forecast errors in the prediction equation.
# The forecasting equation is constructed as follows.  
# First, let y denote the dth difference of Y, which means:

# If d=0:  yt  =  Yt

# If d=1:  yt  =  Yt - Yt-1

# If d=2:  yt  =  (Yt - Yt-1) - (Yt-1 - Yt-2)  =  Yt - 2Yt-1 + Yt-2
#  d=2 is NOT yt = (Yt - Yt-2)
# Note that the second difference of Y (the d=2 case) is not the difference from 2 periods ago.  Rather, it is the first-difference-of-the-first difference, which is the discrete analog of a second derivative, i.e., the local acceleration of the series rather than its local trend.


