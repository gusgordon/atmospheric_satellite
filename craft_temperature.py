"""
Nu = h * L / k

h = Nu * k / L

k = 2.15 for dry air around 240 K

https://www.sfu.ca/~mbahrami/ENSC%20388/Notes/Forced%20Convection.pdf
for Pr > 0.6 laminar:
Nu = 0.664 * Re ** (1 / 2) * Pr ** (1 / 3)
h = (k / L) * 0.664 * Re ** (1 / 2) * Pr ** (1 / 3)

for 60 > Pr > 0.6 turbulent:
Nu = 0.037 * Re ** (4 / 5) * Pr ** (1 / 3)
h = (k / L) * 0.037 * Re ** (4 / 5) * Pr ** (1 / 3)

Re from average Re over wing in each region

Pr = 0.725 for dry air around 240 K

dQ/dt = h * A * deltaT

absortance and emissivity for solar cells ~= 0.8
 -- http://webserver.dmt.upm.es/~isidoro/dat1/Thermooptical.pdf

outgoing:
j* = Stefan_Boltzmann * 0.8 * T ** 4

incoming:
from solar radiation lookup

also from battery & motor heat

https://dothemath.ucsd.edu/2012/01/basking-in-the-sun/

find T of plane

use T to calculate cell efficiency:
http://www.sun-life.com.ua/doc/sunpower%20C60.pdf



positive means heat going into wall

dQ/dt = h * A_wing * (T_air - T_wall) + engine_power * inverse_eff + (solar_irrad * wing_absorbance * A_wing)
        + (albedo * wing_absorbance * A_wing) - (2 * A_wing * Stefan_Boltzmann * 0.8 * T_wall ** 4)
dQ/dt = 0
"""

from scipy import optimize, constants

# When the plane is flying through the air, what will its temperature be?
# After running this and playing around with it, it was obvious that the plane temperature
# is pretty much equivalent to the temperature of air it's flying through
# even accounting for heating due to air compression in front of the wing, heating from waste heat,
# and sun radiation.
def get_craft_temperature(wing_area_lam, wing_area_turb, Re_lam, Re_turb, char_length,
                          T_air, heat_waste, solar_irrad, albedo,
                          absorbance=0.8, emissivity=0.8): # Defaults for wing
    # https://www.sfu.ca/~mbahrami/ENSC%20388/Notes/Forced%20Convection.pdf
    # http://webserver.dmt.upm.es/~isidoro/dat1/Thermooptical.pdf
    # Assume entire aircraft is of same temperature

    initial_guess = 240 # Guess 240 K for initial temp
    air_thermal_conductivity = 2.15 # for dry air around 240k
    Pr = 0.725 # Prandtl number for dry air around 240 K
    
    wing_area = wing_area_lam + wing_area_turb

    # Find heat transfer coeff using Nussult number
    h_lam = (air_thermal_conductivity / char_length) * 0.664 * Re_lam ** (1 / 2) * Pr ** (1 / 3)
    h_turb = (air_thermal_conductivity / char_length) * 0.037 * Re_turb ** (4 / 5) * Pr ** (1 / 3)
    h_eff = (h_lam * wing_area_lam + h_turb * wing_area_turb) / wing_area

    incoming_radiation = solar_irrad * absorbance * wing_area + albedo * absorbance * wing_area

    temp_de = lambda x: (h_eff * wing_area * (T_air - x) + heat_waste + incoming_radiation -
                         (2 * wing_area * constants.Stefan_Boltzmann * emissivity * x ** 4))
    return optimize.fsolve(temp_de, initial_guess)[0]

# Confirm craft temp ~= air temp
assert(239 < get_craft_temperature(7.0, 3.0, 2e5, 7e5, 1.1, 240.0, 50, 1000, 100) < 241)

# This can also be used to get the power requirement for heating battery pack
# Assume insulation (aerogel etc.) blocks all convection and conduction
# Radiation shield (aluminum foil) has emissivity of 0.05
# 400 18650 cells means a pack surface area ~0.35 m^2
# This shows heating takes about 20 W of power
# Due to square cube law, (num_cells / 400) ^ (2 / 3) * 20 = power requirement for heating
#get_craft_temperature(0.0, 0.35, 0, 0, 0.5, 240.0, 20, 0, 0, emissivity=0.05)