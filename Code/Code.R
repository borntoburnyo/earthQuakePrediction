library("dplyr")
library("ggplot2")
library("maptools")
library("astsa")
library("leaflet")


# load dataset 
# dataset download of Kaggle

earthquake <- read.csv("earthQuake.csv", as.is = T, header = T) 

# some clean-ups 

# data-time type
earthquake$Date <- as.POSIXct(strptime(earthquake$Date, "%m/%d/%Y"))

# check NA 
colSums(sapply(earthquake, is.na))
NA.Date <- earthquake[is.na(earthquake$Date), ]

NA.Date$Date <- sapply(NA.Date$Time, function(x) unlist(strsplit(x, "T"))[1])
NA.Date$Time <- sapply(NA.Date$Time, function(x) unlist(strsplit(x, "T"))[2])

NA.Date$Date <- as.POSIXct(NA.Date$Date)
NA.Date$Time <- sapply(NA.Date$Time, function(x) substr(x, 1, 8))

earthquake <- rbind(earthquake[!is.na(earthquake$Date), ], NA.Date)
earthquake <- earthquake[order(earthquake$Date), ]

# attach year/month/day as group variables 
date.earthquake <- sapply(as.character(earthquake$Date), function(x) unlist(strsplit(x, "-")))

earthquake$Year <- as.numeric(date.earthquake[seq(1, 3*nrow(earthquake), by = 3)])
earthquake$Month <- as.numeric(date.earthquake[seq(2, 3*nrow(earthquake), by = 3)])
earthquake$Day <- as.numeric(date.earthquake[seq(3, 3*nrow(earthquake), by = 3)])

# focus on Earthquake Type only 
earthquake <- subset(earthquake, earthquake$Type == "Earthquake")



# Visualiztion 

# Let's look at how many earthquakes happened every year. 
# Overall, the yearly earthquakes has an increasing trend for the past 52 years.    

earthquake %>% 
	group_by(Year) %>% 
  	summarise(Avg.num = n(), 
  		Avg.mag = mean(Magnitude, na.rm = T)) %>%
  			ggplot(aes(x = Year, y = Avg.num)) + 
  				geom_col(fill = "blue") + 
  				stat_smooth(col = "red", method = "loess") + 
  				labs(x = "Year",
  					y = "Total Observations Each Year",
       				title = "Total Observations Each Year (1965-2016)",
       				caption = "Source: Significant Earthquake 1965-2016 by USGS") + 
       			theme_bw()
       			

# Check out the average magnitude of all earthquakes happened each year. 
# The average magnitude tends to be flat around 5.8 approximately, 
# probably this is average over a whole year and extreme strong 
# earthquakes are rare events.   

earthquake %>% 
	group_by(Year) %>% 
		summarise(Avg.num = n(), Avg.mag = mean(Magnitude, na.rm = T)) %>%
			ggplot(aes(x = Year, y = Avg.mag)) + 
  				geom_col(fill = "blue") + 
  				labs(x = "Year",
       				y = "Average Magnitude Each Year",
      				 title = "Average Magnitude Each Year (1965-2016)",
       				 caption = "Source: Significant Earthquake 1965-2016 by USGS") + 
  				theme_bw()
  

# try apply time series models 

# We would frequently prefer to analyze a stationary sequence,
# as this allow us better estimate autocorrelation and other quantities. 
# By saying stationary, we refers to weakly stationarity, where the mean is 
# finite and the same for all $t$, and $\rho(s,t)$ only depends 
# on the lag ($h = \left|s - t \right|$).      

# Now, let's dive in to the time series properties of the data with respect 
# to the average magnitude and total number of earthquakes happened every year. 
# The following time series plot of yearly number of observations shows an overall
# increasing trend, the variance tend to be increasing slightly.
# Let's try differencing first to make it stationary.       

obs <- earthquake %>% 
	group_by(Year) %>%
  	summarise(Avg.obs = n()) 

obs <- ts(obs$Avg.obs, start = 1965, end = 2016)

plot.ts(obs, xlab = "Year", ylab = "Total Observations Yearly", 
	main = "Time Series Plot of Total Observations Yearly")
legend("bottomright", 
	legend = c("Source: Significant Earthquake 1965-2016 by USGS"), 
	cex = 0.6)


# The differenced data seems without trend now, we'll proceed to model part.   

diff.obs <- earthquake %>% 
	group_by(Year) %>%
  	summarise(Avg.obs = n())

diff.obs <- diff(diff.obs$Avg.obs)
diff.obs <- ts(diff.obs, start = 1966, end = 2016) 

plot.ts(diff.obs, xlab = "Year", 
	ylab = "1st Order Difference of Yearly Observations", 
	main = "Differenced Yearly Observation") 
legend("bottomright", 
	legend = c("Source: Significant Earthquake 1965-2016 by USGS"), 
	cex = 0.6)

# After obtain a stationary time series, we could move on to fit possible models. 
# Basic time series models include white noise, moving average (MA), 
# auto regressive (AR), and the combination of the two with/without 
# differencing and seasonality (e.g, ARMA, ARIMA, SARIMA).     

# Usually autocorrealtion and partial autocorrelation plots (ACF/PACF) 
# are the resorts to determine model parameter at first step. 
# ACF measure the correlation of two time points without considering
# the points between the interval, while PACF is autocorrelation measured 
# conditional on all the points between the interval. Here are some properties
# of the basic models, without proof: 1) MA model shows cut-off behavior at
# ACF and tail-off at PACF, 2) AR model shows the opposite behavior. 
# These two are the key points determining model parameters for this study.     

# Let's focus on ACF/PACF plots pick possible model parameters. 
# Apparently, the ACF cut-off at $h = 1$ and the PACF shows cut-off at $h = 2$. 
# We'll start from Auto Regressive Model with differencing ($ARIMA(2,1,1)$).  

acf(diff(obs), main = "ACF of Differenced Data")
pacf(diff(obs), main = "PACF of Differenced Data")


obs.1 <- sarima(obs, 2, 1, 1, details = F)


# Note the $AR2$ paramter is not significant. So, we'll try $ARIMA(1,1,1)$ model. 

obs.2 <- sarima(obs, 1, 1, 1, details = F)

# All the parameters are significant, standardized residuals plot 
# shows no outliers, ACF and LBP test of residuals indicating that 
# the dependencies have been explained by the model, and Q-Q plot shows
#  normal behavior of residuals. It seems $ARIMA(1,1,1)$ model is 
# the *best* one at this moment. The mathematical form of the model is:   

# X_t = 4.25 + 0.42*X_{t-1} + W_t - W_{t-1}

# To have a prediction of the average magnitude of the coming years, 
# we'll follow the same modeling procedure in the previous part. 
# Let's first check the time series plot of the yearly average magnitude. 
# There is a big drop in average magnitude between 1970 and 1980, 
# let's try differencing to make the series stationary. 
# The differenced data looks acceptable, we'll move on to check ACF/PACF.     


mag <- earthquake %>%
	group_by(Year) %>%
  	summarise(Avg.mag = mean(Magnitude))  

mag <- ts(mag$Avg.mag, start = 1965, end = 2016)
diff.mag <- diff(mag)

plot.ts(mag, xlab = "Year", 
	ylab = "Average Magnitude", 
		main = "Time Series Plot of Average Magnitude")
legend("bottomright", 
	legend = c("Source: Significant Earthquake 1965-2016 by USGS"), 
		cex = 0.6)


plot.ts(diff.mag, xlab = "Year", 
	ylab = "Differenced Average Magnitude", 
	main = "Differenced Time Series Plot of Magnitude")
legend("bottomright",
	legend = c("Source: Significant Earthquake 1965-2016 by USGS"), 
	cex = 0.6)


# Both ACF/PACF tails off slowly, it seems $ARIMA(1,1,1)$ would be a good start point.   


acf(diff(mag), main = "ACF of Differenced Data")
pacf(diff(mag), main = "PACF of Differenced Data")


# The diagnostic plots looks good, with white noise behaved residuals, 
# LBP test results show no correlation between residuals, 
# the only point that might need concern is the possible outlier shows 
# in the standardized residuals plot.   


mag.1 <- sarima(mag, 1, 1, 1, details = F)


# Let's also try fit a *bigger* model ($ARIMA(2,1,1)$) to test necessity. 

mag.2 <- sarima(mag, 2, 1, 1, details = F)


# Note, the $AR2$ parameter is not significant. 
# So, the final model we'll use for forecasting is $ARIMA(1,1,1)$, and of the form:  

# X_t & = -0.0026 - 0.86*X_{t-1} + W_t + 0.57*W_{t-1}

# Prediction and conclusion 

# The forecasting of total earthquakes and average magnitude 
# in the next 5 years could be seen in the following plots. 
# There is an increasing trend in total number of earthquakes for 
# the following 5 years. The actual prediction values are 522, 547, 
# 560, 568 and 574. The prediction interval is tight, which shows good 
# confidence. For prediction of average magnitude, the overall trend is
# slightly decreasing, prediction interval looks wide which covers half 
# of the range of magnitudes for the paste 52 years.   

# Typically forecasting of time series models will not be reliable, 
# especially for long term prediction. As the model only using past 
# events as predictors, other predictors are not taken into consideration.
# By considering this, we could conduct analysis using *Vector AutoRegressive Model (VAR)* 
# if we have multivariate time series (other related covariates). 
# It assumes the response values depend not only past value of response time series, 
# but also past values of other timer series.    


par(mfrow = c(2,1))
obs.fore <- sarima.for(obs, 5, 1, 1, 1)
mag.fore <- sarima.for(mag, 5, 1, 1, 1)

