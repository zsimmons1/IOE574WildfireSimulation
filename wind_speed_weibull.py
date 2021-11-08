import pandas
import matplotlib.pyplot as plt
import numpy

data = pandas.read_csv("weather.csv")
wind_speed = data["HourlyWindSpeed"]

#plot of actual data#
N = len(wind_speed)
plt.hist(wind_speed, density=True, alpha=0.5)
plt.title("Hourly Wind Speed over a week (n={})".format(N))
plt.xlabel("Wind speed in m/s")

#Weibull distrbution curve#
from scipy.stats import weibull_min
plt.hist(wind_speed, density=True, alpha=0.5)
shape, loc, scale = weibull_min.fit(wind_speed, floc = 1)
x = numpy.linspace(wind_speed.min(), wind_speed.max(), 100)
plt.plot(x, weibull_min(shape, loc, scale).pdf(x))
plt.title("Weibull fit on win speed data")
plt.xlabel("Wind Speed in m/s")
print(shape, loc, scale)

#generating random wind speeds from a weibull distribution using parameters generated from actual data#
y = weibull_min.rvs(shape, loc, scale, size=N)
plt.hist(y, density=True, alpha=0.5)
x = numpy.linspace(y.min(), y.max(), 100)
plt.plot(x, weibull_min(shape, loc, scale).pdf(x))
plt.title("Hourly Wind Speed over a week (n={})".format(N))
plt.xlabel("Wind speed in m/s")
