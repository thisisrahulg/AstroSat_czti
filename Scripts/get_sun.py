from astropy.time import Time
from astropy.coordinates import get_sun


t = input(" Enter time in UTC [ Y-m-dTh:m:s ] : ")

time = Time(t, format = 'isot', scale = 'utc')

print ( get_sun(time))

