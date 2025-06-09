"""
unit conversion
Note that all the unit in computation is in SI unit system, unless specified.
"""

# speed conversions
const g_kmh_per_ms = 3.6
const g_mph_per_ms = 2.23694
const g_seconds_per_hour = 3600.0
const g_j_per_kwh = 3.6e6
const g_diesel_j_per_liter = 38.6e6 #https://en.wikipedia.org/wiki/Energy_density
const g_liter_per_gallon = 3.78541
const g_diesel_j_per_gallon = g_diesel_j_per_liter * g_liter_per_gallon
const g_meter_per_mile = 1609.344

global g_longhaul_distance = 500.0 * g_meter_per_mile

ms2kmh(spd::Real) = (spd*g_kmh_per_ms)
ms2mph(spd::Real) = (spd*g_mph_per_ms)
kmh2ms(spd::Real) = (spd/g_kmh_per_ms)
mph2ms(spd::Real) = (spd/g_mph_per_ms)
second2hr(s::Real) = (s/g_seconds_per_hour)
hr2second(hr::Real) = (hr*g_seconds_per_hour)
j2kwh(j::Real) = (j/g_j_per_kwh)
kwh2j(kwh::Real) = (kwh*g_j_per_kwh)
j2liter(j::Real) = (j/g_diesel_j_per_liter)
liter2j(l::Real) = (l*g_diesel_j_per_liter)


