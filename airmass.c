/*******************************************************************************
 *
 *  airmass -- A program to compute the airmass for a given zenith angle.
 *
 *  Copyright (C) 2000 Reed D. Meyer (rxm128@yahoo.com).
 *  You are permitted to copy and distribute this program as much as you want,
 *  as long as you leave the source code, including this comments section,
 *  intact and unmodified.
 *
 *  To compile, type:  gcc airmass.c -o airmass -s -O3 -lm
 *     (This command should work in any environment where the GNU C compiler,
 *     gcc, is installed, including Linux and most UNIX environments)
 *
 *  To run, type:
 *
 *    airmass [-a <alt>] [-d <day>] [-H <h>] [-l <lat>] [-L <wavelength>]
 *            [-p <P>] [-t <T>] [-v] <zenith>
 *
 *  where:
 *     <zenith> is the APPARENT zenith angle (not the true angle in the
 *         absence of atmosphere) at which you want the airmass;
 *     -a <alt> sets the altitude of the observer to <alt> meters above mean
 *         sea level (the default is 0 meters);
 *     -d <day> specifies the date, in days since the beginning of the year
 *         (0 = midnight Jan. 1; the default is 80, corresponding roughly to
 *         the vernal equinox); don't worry about fractions of a day, since the
 *         seasonal effect on airmass is relatively small, as discussed below;
 *     -H <h>, if specified, tells the program to additionally compute an
 *         estimate of the water vapor column density and airmass, using a
 *         CRUDE model for the water vapor profile based on the local relative
 *         humidity at ground level, <h>.  <h> is in percent (100 = saturated).
 *         The local humidity is NOT used in the computation of the standard
 *         (dry atmosphere) airmass.  See EXTINCTION below for details.
 *     -l <lat> sets the observer's latitude to <lat> degrees (negative = south;
 *         default is 45 degrees);
 *     -L <wavelength> changes the observing wavelength in Angstroms (default
 *         5500 A; don't worry too much about an exact figure here, as the
 *         wavelength dependence is very small, as discussed below);
 *     -p <P> specifies the local pressure in millibars, presumably read off
 *         a barometer at the time of observation (default is 1013.25 mb,
 *         the so-called standard atmosphere);
 *     -T <T> specifies the local temperature in Celsius, again at the time
 *         of observation (default is 15 C, again the standard atmosphere);
 *     -v is a flag that increase the "verbosity" of the output.  In this
 *         context, the program computes the relative improvement it makes
 *         over some standard airmass formulas astronomers use (e.g. polynomial
 *         expansion in (sec(z) - 1) terms) -- see VERBOSITY OPTION below
 *         for details.
 *
 *
 *  VERSION HISTORY
 *  ===============
 *
 *  1.0 (2000 Aug. 13): first public release.
 *
 *
 *  BACKGROUND
 *  ==========
 *
 *  This program, and the necessary underlying algorithm, was developed to
 *  find out to just what degree the standard approximate formulas for airmass
 *  are inaccurate.  Many astronomers use the following polynomial in
 *  (sec(z) - 1) to calculate the airmass for a given zenith angle z:
 *
 *                                                     2                      3
 *  X = secz - 0.0018167(secz - 1) - 0.002875(secz - 1)  - 0.0008083(secz - 1)
 *
 *  (that is, those astronomers who do not use the simpler approximation
 *  X = sec(z), or the formula used in IRAF which is from Allen's _Astrophysical
 *  Quantities_.)  The fact is, though, that some astronomers take this formula
 *  practically as "truth" and do not stop to think about its origin.
 *  This formula is a fit to a table of theoretical airmass as a function of z,
 *  computed at the turn of the century, when the structure of the atmosphere
 *  was less well-understood and equations were simplified to ease numerical
 *  integration since there wasn't the benefit of computers.  To the best of
 *  my knowledge, the polynomial approximation above was first presented by
 *  R.H. Hardie in _Astronomical Techniques_ (W.A. Hiltner, ed.) (Chicago:
 *  U. Chicago Press), p. 180 (1962).  This was a fit to computations by
 *  A. Bemporad (Mitt. d. Grossh. Sternwarte Heidelberg 4 (1904)), published
 *  again by E. Schoenberg in _Handbuch der Astrophysik_ (Berlin: J. Springer),
 *  2, 268-273 (1929).*
 *       So there are two potential sources of inaccuracy in using this method:
 *  First, and perhaps the greatest source, is simply that the formula is
 *  a fitted function to a table of discrete values (and, moreover, was not
 *  designed to fit well the entire range of zenith angles).  The second is
 *  that the underlying table of theoretical values itself is based on an
 *  old model atmosphere, and used simplifying assumptions which further
 *  degrade the accuracy of the model.
 *       What if, instead, we were to use the power of modern computers to
 *  obviate the need for most simplifying assumptions and, in addition,
 *  calculate the airmass directly for any given z rather than make use of
 *  interpolating or best-fit functions?  This is what this program was
 *  designed to do.  In carrying out such a "modernization" of the airmass
 *  function, we might as well incorporate more modern theories of the
 *  structure of the atmosphere, too.
 *
 *  *Interestingly, the literature seems confused by Bemporad's work: some
 *   think the (sec(z) - 1) fit was to real DATA but in fact Bemporad just
 *   tabulated the results of his theoretical airmass formula.  It is pleasing
 *   to note that even the old German references can be found in the Yale
 *   Astronomy Library.
 *
 *
 *  DERIVATION
 *  ==========
 *
 *       We assume a spherically symmetric atmosphere of constant chemical
 *  composition.  The equation for the column density, which can be derived
 *  strictly from Snell's Law with the above assumption, is (c.f.
 *  Schoenberg (1929)):
 *
 *                 /\ r_m
 *                 |                rho dr
 *                 |     -----------------------------
 *        N(z) =   |                 2     2      2                       (1)
 *                 |                r_o  mu_o  sin  z
 *                 |     sqrt(1 -   -----------------)
 *                 |                    2   2
 *                \/ r_o               r  mu
 *
 *  where the subscript o refers to values at the observer, and r is the
 *  distance from the center of the earth, mu is the air's index of refraction
 *  at r, rho is the air's density at r, r_m is the upper limit on the
 *  integration (taken as the top of the mesosphere; see below) and z
 *  is the apparent zenith angle for which we're calculating the airmass.
 *  Note that this is just the expression for the column density of air
 *  directly overhead, with a correction factor (sqrt(...) in the denominator)
 *  to account for the light's actual path through the atmosphere when the
 *  source is not actually at the zenith.
 *       To get the airmass as it's traditionally defined (i.e., normalized to
 *  1 at the zenith), simply divide by the column density at z = 0:
 *
 *                   N(z)
 *           X(z) = ------  .                                             (2)
 *                   N(0)
 *
 *  The Clausius-Mossotti equation (see, for example, B. Garfinkel, AJ 72, 235
 *  (1967)) gives the index of refraction in terms of the density to an
 *  extremely good approximation; rewriting the equation gives
 *
 *             2     3 + 4 c rho
 *           mu  =  -------------                                         (3)
 *                   3 - 2 c rho
 *
 *  where c is a wavelength-varying constant depending only on the chemical
 *  composition of the gas.  A formula for air's index of refraction at a
 *  given wavelength appearing in the CRC Handbook of Chemistry and Physics
 *  can be rewritten to yield c:
 *
 *                      13.412      0.3777      -7     273.15 K
 *      c ~ [2875.66 + -------- + ----------] 10   R [-----------]        (4)
 *                      2   -8      4   -16            101325 Pa
 *                     L  10       L  10
 *
 *  where L is the wavelength in Angstroms and R is the ideal gas constant.
 *  (c has units of volume/mass so R needs to be expressed in energy/mass/K).
 *
 *  The main problem is finding an expression for rho(r).  We do this by
 *  assuming a temperature profile and then combining the ideal gas law
 *  (which expresses density in terms of pressure and temperature) with the
 *  condition of hydrostatic equilibrium (which gives density in terms of
 *  pressure).  This yields a unique density profile, given an assumed
 *  temperature profile.  The problem is the temperature profile.  This can
 *  only be determined uniquely by solving the equations of fluid dynamics
 *  for the atmosphere coupled with knowledge of the radiative processes in
 *  the atmosphere with appropriate boundary conditions, i.e., the incredibly
 *  complex problem which still makes accurate weather forecasting impossible.
 *  So we assume a temperature profile which is based on average measured
 *  values of temperatures at different heights in the atmosphere.  As I will
 *  argue below under ASSUMPTIONS, I feel that the approximate temperature
 *  profile probably does not deviate strongly from the real profile at the
 *  observer under most conditions, and since this is the primary assumption
 *  made in calculating the airmass for dry air, said airmass probably doesn't
 *  differ strongly from the real value.
 *       For our temperature profile, we start by assuming the atmosphere is
 *  broken down into several layers within which the temperature changes
 *  linearly with altitude.  In any given layer n, the temperature is
 *  T(Q) = T_n + beta_n (Q - Q_n), where beta_n is the slope of the temperature
 *  rise with altitude within the layer n, T_n is the temperature at the bottom
 *  of that layer, Q is the height (defined below), and Q_n is the height of
 *  the bottom of the layer.  This is the methodology used in the _U.S. Standard
 *  Atmosphere, 1962_ (Washington, DC: U.S. Govt. Printing Office), wherein a
 *  "standard atmosphere" appropriate for the continental United States is
 *  defined (a "standard atmosphere" being useful for aircraft and so forth).
 *  We take their bottom eight temperature layers (corresponding roughly to
 *  the troposphere, stratosphere, and mesosphere; top of eighth layer is at
 *  90 km above mean sea level).  But in order to generalize their work to make
 *  it more appropriate for any latitude and any season, we change the values
 *  defining the lowest two layers (the troposphere and lower stratosphere).
 *  We let the boundary between those two layers (the tropopause; Q_1) vary as
 *  a function of time of year and latitude.  The same thing goes for beta_0,
 *  the temperature falloff in the first layer.  The bottom of the first layer
 *  (Q_0) is set to the height Q of the observer's location.  The boundary
 *  condition "temperature at bottom of first layer" (T_0) is set to the
 *  actual local temperature at the observer.  The condition T_1 equals the
 *  temperature at the top of the first layer, i.e., as determined from the
 *  above expression for T with n = 0 and Q = Q_1.  Finally, beta_1 is set
 *  by forcing the temperature at the top of the second layer to match T_2
 *  (fixed at 216.65 K, the value in the U.S. Standard Atmosphere).  Although
 *  the temperatures in the third layer and higher (altitudes > 20 km) *do*
 *  vary with season and latitude, the variations are less severe than in the
 *  troposphere and lower stratosphere, and fortunately for us the atmosphere
 *  has already greatly thinned out by that point so these layers contribute
 *  little to the airmass calculation anyway.
 *       To determine Q_1 and beta_0 (the two quantities allowed to vary with
 *  time and latitude, and which define the lower two layers along with T_0),
 *  I extracted a table of values from two plots appearing on pages 48-49 of
 *  _General Meteorology_, 3rd ed.,  by H.R. Byers (New York: McGraw-Hill)
 *  (1959).  Byers provides a plot of the temperature distribution within the
 *  troposphere as a function of latitude for summer, and another for winter.
 *  Keep in mind that these are average values.  (In particular, the location
 *  of the tropopause at any given time is far more complex than the smooth
 *  line drawn in Byers' plots, due to the movement of the jet streams.)  I
 *  then fit a polynomial to the data to get a smooth function for beta_0 for
 *  any given latitude.  I fit another polynomial to the data to get a function
 *  for "r1", the altitude of the troposphere (which is related to Q_1).  The
 *  polynomials fit quite well and in any case are far more than adequate given
 *  the fact that Byers' plots are time averages.  beta_0 needed only a
 *  4th-order polynomial to fit the data but r1 required a 10th-order
 *  polynomial.  You can find the values of the best-fit coefficients by
 *  looking at the declarations of betapoly[] and r1poly[] in the source code.
 *  (The units are K/m and km, respectively.)  Then, to simulate the effect of
 *  seasonal changes, the odd coefficients in the polynomials sinusoidally vary
 *  with day of the year (period = 1 year; no change at day number 202,
 *  corresponding to the peak of summer).
 *       Now that we have values for the beta_n's, T_n's, and Q_n's (either
 *  from the _U.S. Standard Atmosphere_ or the method above for the bottom two
 *  layers), the simplicity of the temperature model leads to easily solvable
 *  equations which result in a closed-form expression for the density:
 *
 *            /
 *           |            T_n            (1 + g_E /(R beta_n))
 *           | [------------------------]
 *   rho     |   T_n + beta_n (Q - Q_n)                      , beta_n <> 0
 *  -----  =_/                                                              (5)
 *  rho_n    \
 *           |        g_E
 *           | exp[ ------- (Q_n - Q) ]
 *           |       R T_n                                   , beta_n = 0
 *            \
 *
 *  where as before, subscripts n refer to the bottom of some temperature layer
 *  n (0 <= n < 8); g_E is the standard value of the acceleration of gravity
 *  at the earth's surface (9.80665 m/s^2); and all other symbols are as
 *  defined above.
 *       In order to use (5) in Equation (1), we need to write r in terms of Q.
 *  (The only other quantity in the integrand of (1) which varies with r is mu,
 *  and we already have an expression for that in terms of rho.)  So, let us
 *  discuss the definition of Q.  The _U.S. Standard Atmosphere_ defines the
 *  geopotential altitude (here referred to as "h") as the integral from 0 to r
 *  of (g/g_E) dr (compare to their equation I.2.5-(1)), where g is the gravity
 *  at some r.  They then define their eight lowest temperature layers in terms
 *  of the geopotential altitude h.  The advantage of changing variables to h
 *  is that it lets us forget the fact that the g field is varying, thus
 *  simplifying the calculations.  To write the equations more simply, I
 *  further define the quantity Q by Q = h - (r_E^2 / r_msl), where r_E is the
 *  "standard" radius of the earth (i.e., the radius you would get from the
 *  formula g_E = G M_E / r_E^2, where M_E is the mass of the Earth, G is the
 *  gravitational constant, and g_E is the "standard" value of gravity as
 *  before); and r_msl is the actual radius of the earth at the latitude of
 *  the observer (i.e., the distance between the earth's center and a point
 *  at mean sea level at the latitude of the observer; "msl" means "mean sea
 *  level").  r_msl varies because the earth is an oblate spheroid; the reason
 *  for distinguishing r_E and r_msl is that Earth's gravity at the equator is
 *  slightly less than at the poles, because the earth's radius is larger at
 *  the equator.  Note that Q is always negative.
 *       We ignore the fact that the gravity vector is not exactly equal to
 *  the gravity vector of a spherical, nonrotating Earth (see ASSUMPTIONS
 *  below).  Combining our definitions for h and Q lets us write r in terms
 *  of Q:
 *
 *                      2
 *                   r_E
 *            r = - -----                                                 (6)
 *                    Q
 *
 *  Similarly we can find a Q for any desired r; for example, Q_o (the lower
 *  limit on the integrand in Equation (8), and the Q of the observer) is
 *  related to the observer's altitude "a" above mean sea level by
 *  Q_o = -r_E^2 / (r_msl + a).
 *       We calculate r_msl from the analytic geometry of the oblate spheroid,
 *  and the standard definition of the geographic latitude.  (Latitudes of
 *  actual places on the earth are always geographic latitudes.)  The result is
 *
 *              4     4    4     2
 *      2      b  + (a  - b ) cos  phi
 *     r    = -------------------------                                   (7)
 *      msl     2     2    2     2
 *             b  + (a  - b ) cos  phi
 *
 *  where a is the major axis of the spheroid (the Earth's equatorial radius),
 *  b is the spheroid's minor axis (the Earth's polar radius), and phi is the
 *  observer's latitude.  (Incidentally, the geocentric latitude phi_c can be
 *  found from:  tan phi_c = (b^2 / a^2) tan phi.)
 *       Plugging (6) into (1) gives
 *
 *                 /\ Q_m                  2
 *                 |                rho r_E  dQ
 *                 |     -------------------------------
 *        N(z) =   |                    2   2      2                      (8)
 *                 |      2            Q  mu_o  sin  z
 *                 |     R  sqrt(1 -  -----------------)
 *                 |                       2     2
 *                \/ Q_o                  Q_o  mu
 *
 *  (Q_m being the Q at the top of the highest temperature layer; for this
 *  layer, h = 88.743 km.  We take this as a suitable upper limit on the
 *  integration because the mass above this height is negligible; the air
 *  density at this altitude is 5 orders of magnitude less than at sea level.
 *  A practical reason for choosing this height for our upper limit is that
 *  at about this height (the top of the mesosphere), the molecular weight of
 *  the air begins to change appreciably.  While we can safely assume that all
 *  the air below this altitude has a fixed proportion of all the components
 *  important for extinction at visible and near-infrared wavelengths (besides
 *  water vapor, ozone, and aerosols), this assumption is not a good one
 *  above this rough altitude.
 *       Plugging Equations (3), (4), and (5) into (8) and integrating gives
 *  N(z); then we set z = 0 to use the same formulas and integration routines
 *  to get N(0); Equation (2) then gives us the desired airmass at z.
 *
 *
 *  EXTINCTION
 *  ==========
 *
 *  Of course, the whole point of calculating airmass, for astronomical
 *  purposes anyway, is to determine the atmospheric extinction.  In addition
 *  to calculating the airmass, the program also computes the column density,
 *  because it is perhaps a more appropriate figure for determining extinction.
 *  After all, the airmass only tells you what the extinction is going to be
 *  at some zenith angle, RELATIVE to at the zenith; you still need to know
 *  the extinction at the zenith!  The standard approximate formulas for
 *  airmass cannot provide this important bit of information, which is directly
 *  related to the column density.
 *       Several sources contribute to the atmospheric extinction at any given
 *  wavelength.  For visible and infrared wavelengths, we can limit the
 *  discussion to four primary contributors; the others are small enough to
 *  be negligible.  (Actually, there is an EXTREMELY important fifth component,
 *  namely suspended water droplets, commonly known as clouds; but in this
 *  discussion we'll assume that you're not trying to do high-precision
 *  photometric or spectroscopic work through clouds, so we'll ignore them.)
 *  The four components are: dry air, water vapor, ozone, and dust.  "Dust"
 *  means any particulate matter suspended in the air (i.e., aerosols), and
 *  "dry air" signifies clean (dust-free) air with zero humidity, and having
 *  a constant chemical composition.*
 *       The total extinction, in magnitudes, is given approximately by:
 *
 *      E = 2.5 log(e) (k_air N_air + k_O3 N_O3 + k_dust N_dust + k_H20 N_H20)
 *
 *  where the various k's are the extinction coefficients, in cm^2/g, for each
 *  contribution; the N's are the column densities (in g/cm^2) of each
 *  contribution; and the subscripts denote the contributions: air = dry,
 *  clean air of constant (or negligible) ozone composition; O3 = ozone;
 *  dust = suspended particles of any sort; H20 = water vapor.  The k's are
 *  wavelength-dependent, making E a function of wavelength.  This program
 *  computes the column density N_air, and optionally (with the "-H" flag on
 *  the command line) a crude estimate for N_H20.  (If your k's are in units of
 *  cm^2 rather than cm^2/g, convert the column densities to number densities
 *  by multiplying by A/m:
 *
 *       N (number density, units 1/cm^2) = N (g/cm^2) * A/m
 *
 *  where A is Avogadro's constant and m is the mean molecular weight of the
 *  constituent.  For dry air, m = 28.9644 g/mol, so the number density equals
 *  2.07915 * 10^22 * the column density in g/cm^2 for the dry air component.)
 *  k_air * N_air is approximately 0.2 at visual wavelengths, for dry air,
 *  observing at the zenith.
 *       Let us discuss the "-H" option briefly and how a crude estimate of
 *  the water vapor airmass and column density is calculated.  It is impossible
 *  to calculate the water vapor values from a simple theory; it requires the
 *  same sort of full-blown analysis that makes the calculation of the
 *  temperature profile impossible.  Thus, one needs either a direct
 *  measurement of the water vapor as a function of altitude (for example,
 *  through LIDAR), or one must make some approximation.  Based on an
 *  admittedly brief look into real vapor profiles, I propose the following
 *  extremely simplistic model, which does have the desired feature of being
 *  constrained by the actually measured value of the humidity at the
 *  observatory:
 *
 *            /
 *           |  H_0                                          , Q_0 < Q < Q_1
 *           |
 *           |          Q_2 - Q
 *     H =  _/  H_0 [ ----------- ]
 *           \         Q_2 - Q_1                             , Q_1 < Q < Q_2
 *           |
 *           |
 *           |  0                                            , otherwise
 *            \
 *
 *  where the Q's are as defined in the section "DERIVATION", H is the relative
 *  humidity at some point Q, and the H_n's are the relative humidities at the
 *  Q_n's.
 *       The primary motivation behind the above linear model (constant H up
 *  to the troposphere; linear falloff through the lower stratosphere; zero
 *  humidity above that) is that, as I just said, it is constrained by the
 *  measured value of the humidity at ground level (H_0).  In fact, it is
 *  quite strongly constrained, since H throughout the troposphere is assumed
 *  equal to this value.  Believe it or not, even though this is a very lame
 *  model, I believe it does have some merit, for the following reason: looking
 *  at a couple typical vapor profiles, the relative humidity did not tend to
 *  vary THAT MUCH throughout the lower troposphere -- "THAT MUCH" meaning
 *  "deviations larger than a factor of two, over large regions of the
 *  troposphere, seemed rare".  Furthermore, there is one great numerical
 *  property of the relative humidity which "limits the carnage" as far as
 *  keeping the REAL relative humidity from really deviating by orders of
 *  magnitude from the assumed value of H_0.  Which is this: a relative humidity
 *  above 100% is forbidden.  If the relative humidity exceeds 100%, it will
 *  not be by very much, and will probably not last for long, because the
 *  water vapor will try to condense as soon as it can (e.g., it spots a dust
 *  particle floating by and grabs onto it).  So, for all intents and purposes,
 *  we are numerically limited to values between 0% and 100% everywhere in the
 *  atmosphere, which is still a huge range of possible humidities but at least
 *  it is not all the real numbers.
 *       I am interested in seeing how well the above crude formula for H
 *  holds up against real data (both real measured humidity profiles as well
 *  as water vapor extinction computed from accurate photometry at some
 *  wavelength where water vapor dominates).  If you can supply any feedback
 *  on this, or if you know of a better humidity model which doesn't depend
 *  on complicated fluid dynamics equations, please do let me know.
 *       Anyway, once we have an assumed humidity profile, it is straightforward
 *  to calculate the water vapor column density and airmass.  The definition
 *  of relative humidity gives the water vapor pressure, given the vapor
 *  pressure at saturation.  The Clausius-Clapeyron equation gives the
 *  saturation vapor pressure as a function of temperature.  For the Clausius-
 *  Clapeyron equation, see Equation 4.4 in A.S. Monin, _Weather Forecasting
 *  as a Problem in Physics_ (Cambridge, MA: MIT Press), p. 15 (1972); we
 *  rewrite Monin's formulation of the equation so that it gives us exactly
 *  what we want, the saturation vapor pressure as a function of temperature,
 *  relative to some saturation pressure at some standard temperature.  These
 *  two conditions (the humidity profile and the Clausius-Clapeyron equation)
 *  let us solve for the vapor pressure.  Then, we assume local thermodynamic
 *  equilibrium (i.e., that the local temperature of the water vapor equals the
 *  temperature of dry air at that point -- a good assumption) and that the
 *  water vapor pressure behaves independently of the other gases in the
 *  neighborhood (not completely true, but certainly a good enough assumption
 *  given the poorly-known humidity profile).  Finally, we need only apply the
 *  ideal gas law (very close to valid for water vapor) and we have our density
 *  as a function of location in the atmosphere.  We use the same algorithm
 *  as with dry air, from that point on, to compute the column density and
 *  hence the airmass of the water vapor at some zenith angle z.
 *
 *  *In our airmass and column density calculations, ozone is included in the
 *   dry air component (with a constant fraction of the total air pressure
 *   assumed), showing up as a small contribution to the mean molecular weight.
 *   The mean molecular weight is the only way ozone has any impact in the
 *   dry air calculation.  Because the ozone concentration does change
 *   significantly throughout the atmosphere (it peaks in the lower
 *   stratosphere, and anyone who follows the news knows that there are "ozone
 *   holes", i.e., ozone is also a function of latitude and longitude), its
 *   column density CANNOT be calculated with the current code.  You should
 *   determine the ozone extinction contribution from a separate algorithm
 *   and then apply that contribution to the total extinction using the formula
 *   for E above.  (Since ozone is a small contribution to the mean molecular
 *   weight, do not worry about its effect in the calculation of the dry
 *   air column density or airmass.)
 *
 *
 *  NOTE ON ASSUMPTIONS
 *  ===================
 *
 *  Here is a rough analysis of the expected error in the calculation of the
 *  column density and airmass from using various assumptions.  As far as I
 *  know, the following list contains all the assumptions this program makes
 *  use of.  The list is sorted in order of their expected impact for most
 *  zenith angles (the worst assumption listed first).
 *       Note that this assumption list is only for the DRY AIR column density
 *  and airmass.  It is not a list of the assumptions on the water vapor column
 *  density or airmass which the program optionally calculates (see EXTINCTION
 *  above for more on that).  Nor does it deal with ozone or aerosol column
 *  densities or airmasses.
 *       Since the exact quantitative impacts of the assumptions depend on
 *  multiple factors, most particularly the zenith angle of interest, I only
 *  treat them in a qualitative sense, using adjectives like "moderate",
 *  "small", etc. to denote their expected relative importance.  From most to
 *  least in terms of impact, the adjectives are "moderate", "small", "very
 *  small", "minimal", "insignificant", and "unimportant".  Based on some
 *  VERY crude and quick calculations, "moderate" would correspond to a percent
 *  error (from using the approximation) of roughly 1% or less; "small" would
 *  correspond to something like 0.2% or less; "very small" smaller still;
 *  "minimal" meaning something like 0.01% if not considerably less;
 *  "insignificant" here meaning "probably falls into the class 'unimportant',
 *  but I'm not sure, and in any case is probably not a bigger assumption than
 *  something in the 'minimal' class"; and "unimportant" meaning around 0.001%
 *  if not considerably less.
 *       One thing to keep in mind is that these descriptions really only apply
 *  to LOW TO MODERATE ZENITH ANGLES (i.e., z less than ~ 60 degrees).  Many
 *  of the assumptions tend to become worse with z; those that are believed to
 *  have a significant z dependence (so that the assumptions are likely to be
 *  notably worse at extreme zenith angles) have the remark "z dependence."
 *  For example, the temperature profile assumption might become quite severe
 *  (i.e., larger than the rough 1% upper limit corresponding to the so-called
 *  "moderate" class) for zenith angles near 90 degrees.
 *       Another thing to keep in mind is that under no circumstances are the
 *  computed column density or airmass of higher accuracy than about 1 part in
 *  10^8.  This is because we consider the numerical integration "converged"
 *  when the column density changes by less than 1 part in 10^8 between
 *  iterations.  One could always change this tolerance limit, but there seems
 *  little point when certain assumptions like the temperature profile likely
 *  contribute a far larger fractional error than 1 in 10^8.  (The convergence
 *  criterion is the motivation, by the way, for printing the results to 8
 *  significant digits.  If all the following assumptions were negligible, in
 *  principle all 8 digits would be valid.)
 *
 *  TEMPERATURE PROFILE:  SMALL to MODERATE for column density, SMALL for
 *       airmass.  z dependence (possibly strong).  This is almost certainly
 *       the biggest assumption for the dry air calculations.  The problem is
 *       that we cannot compute this profile theoretically, and it would be
 *       impractical for most observatories to launch radiosondes or whatever
 *       to probe the atmosphere overhead.  Fortunately, Nature is kind to us
 *       in that the day-to-day fluctuations in the temperature at a given
 *       altitude are generally small compared to the absolute temperature;
 *       that is, the percent error of the fluctuations is relatively small
 *       when the temperature is expressed in Kelvins.  According to
 *       Figure II.2.1(b) in the _U.S. Standard Atmosphere_, at any given
 *       altitude the temperature will tend to be within about 10 K of the
 *       average value, over much of the altitude range of interest.  Thus
 *       the typical percent error is perhaps 5% at any given altitude and
 *       time under most conditions.  The variance is greatest at ground level,
 *       but we correct for that to a great extent by *incorporating* the
 *       ground temperature explicitly into our calculations of the temperature
 *       profile for the bottom two layers.  This should eliminate much of the
 *       variance at the lowest altitudes.  Furthermore, systematic offsets
 *       in one direction at one range of altitudes get balanced to some
 *       extent by offsets at another range, so the overall effect is lessened.
 *       And, unlike the U.S. Standard Atmosphere, we incorporate seasonal and
 *       latitude trends into our calculation of the bottom two layers, which
 *       helps reduce the offset between the theoretical profile and the real
 *       one.  Based on a crude experiment, I expect the error from not knowing
 *       the true profile to affect the column density by < 0.8%, and the
 *       airmass by < 0.3%.  These are conservative estimates, and the real
 *       error might be considerably less.
 *  SPHERICAL SYMMETRY:  SMALL under most conditions.  z dependence (possibly
 *       strong).  Except for the calculation of r_msl (see the discussion of
 *       the next assumption), complete spherical symmetry is assumed.  This,
 *       of course, is never absolutely true, for otherwise there would not
 *       be, for example, high-pressure and low-pressure weather systems since
 *       all locations on earth would have the same weather conditions at a
 *       given altitude.  The question really is to what extent these deviations
 *       from spherical symmetry affect the true airmass.  These deviations
 *       do not affect the column density at z=0 (the zenith) but in principle
 *       would for any other z.  Weather systems are essentially confined to
 *       the troposphere and lower stratosphere, i.e., roughly the lowest 20 km.
 *       Thus, except for extreme zenith angles, large weather patterns would
 *       not really affect our airmass that much because the horizontal distance
 *       starlight would travel through the weather system would not be great
 *       (less than 20 km for z = 45 degrees).  So, only local deviations from
 *       symmetry are important.  Small-scale effects such as turbulence would
 *       average out over the light's path, so only mesoscale deviations could
 *       be important, i.e., regions with sizes on the order of 20 km * sec(z)
 *       or less, but larger than ~ 1 km.  Furthermore, the dominant impact of
 *       these varying mesoscale regions would be a varying temperature, but
 *       this should be on the same order of magnitude if not considerably
 *       smaller than the vertical temperature variations, and therefore the
 *       effect of this assumption would be about the same if not much smaller
 *       than the "temperature profile" assumption.  One special case of
 *       deviation from spherical symmetry is discussed under the next
 *       assumption heading.
 *  GRAVITY VECTOR:  SMALL.  As the _U.S. Standard Atmosphere_ points out, the
 *       true gravity vector at any point does not generally point straight
 *       down to the center of the earth nor have the magnitude that you would
 *       calculate from Newton's Law of Gravitation for a sphere.  There are
 *       two main reasons.  The somewhat more important one, numerically, is the
 *       centrifugal force caused by the rotation of the earth.  Additionally,
 *       the earth is really an oblate spheroid (ironically, the earth is
 *       oblate precisely BECAUSE of the centrifugal force).  In principle, we
 *       could account for these effects in the computer code, since the theory
 *       is quite well known.  However, there really is no reason to complicate
 *       the algorithm in this manner, because the effects are rather small.
 *       Compared to the formula for a gravitating sphere, the first term of
 *       the centrifugal correction is ~ 580 times smaller, and the first term
 *       of the oblate spheroid correction is ~ 920 times smaller.  These are
 *       both on the order of the error of the "temperature profile" assumption
 *       above.  Thus, it did not seem necessary to complicate the code unless
 *       we have better information on the temperature profile.  The oblateness
 *       of the earth *is* accounted for in one way, however, since it is
 *       easy to implement it in the code.  It is in the calculation of r_msl,
 *       the true distance between mean sea level on some point on the earth's
 *       surface and the earth's center.  (See Equation (7) under DERIVATION
 *       above.)
 *  CONSTANCY OF CHEMICAL COMPOSITION:  VERY SMALL.  z dependence.  This is the
 *       assumption that all the components of dry, clean air have the same
 *       partial pressures throughout the atmosphere region over which we
 *       integrate (i.e., throughout the lowest 90 km).  For all the atmospheric
 *       gases which are important to extinction in visible and near-infrared
 *       wavelengths, except water vapor and ozone, this is a very good
 *       approximation.  It is such a good approximation that the U.S. Standard
 *       Atmosphere assumes a constant molecular weight of 28.9644 throughout
 *       the lowest 90 km.  Water vapor is considered separately (it is not
 *       a component of "dry, clean air", of course); ozone is a special case
 *       in that it is included into the dry air calculations through its
 *       molecular weight, but is considered constant throughout the atmosphere
 *       even though it's not.  (Ozone's contribution to the total molecular
 *       weight is small in any case, no matter how we consider its contribution
 *       to dry air.)  See also the "airmass independent of contaminants"
 *       assumption below.
 *  IDEAL GAS:  MINIMAL.  In many sources in the literature, dry air is said
 *       to behave very similarly to a ideal gas.  The deviation from an ideal
 *       gas is without question smaller than some of the assumptions above.
 *  HYDROSTATIC EQUILIBRIUM:  MINIMAL.  For hydrostatic equilibrium to be
 *       seriously violated, a parcel of air would practically have to be in
 *       a state of free fall.  Small-scale regions of turbulence would average
 *       out over a light's path, so a large region would have to be undergoing
 *       this type of motion extreme.  This seems quite unlikely except maybe
 *       near a tornado.  In which case, you wouldn't be observing.  Thus,
 *       departure from hydrostatic equilibrium is expected to have a quite
 *       minor effect on the true airmass.
 *  AIRMASS INDEPENDENT OF CONTAMINANTS:  INSIGNIFICANT.  z dependence.
 *       "Contaminants" here means water vapor, ozone, or aerosols (where
 *       ozone abundance varies with height).  As stated in the section
 *       "IMPORTANCE OF ATMOSPHERIC VARIABLES", water vapor is a very weak
 *       contribution to the total airMASS.  The same goes for ozone and
 *       aerosols.  We treat those components separately when considering
 *       EXTINCTION (q.v.), since water vapor, ozone, and aerosols do contribute
 *       heavily to the total extinction, even though they contribute so little
 *       to the mass of dry air.
 *  CLAUSIUS-MOSSOTTI EQUATION:  INSIGNIFICANT.  z dependence.  As it is, the
 *       Clausius-Mossotti Equation describes pretty well the relation between
 *       index of refraction and density for many materials, including air (c.f.
 *       Garfinkel 1967).  Furthermore, the index of refraction for air anywhere
 *       within the earth's atmosphere is going to be pretty close to 1, no
 *       matter what the functional form of mu(rho) is.  Finally, mu has only
 *       a rather weak impact on the airmass calculation (mainly because mu
 *       does stay close to 1); it is more important in calculating refraction.
 *       All these facts combined together illustrate that any errors which
 *       arise from using the Clausius-Mossotti Equation to calculate mu are
 *       probably quite small.
 *  ALL MASS BELOW 90 KM:  UNIMPORTANT.  The density at 90 km is about 10^5
 *       times less than at sea level, and continues to fall off nearly
 *       exponentially above that (see e.g. Figure I.2.11 in the _U.S. Standard
 *       Atmosphere_).  Thus, the contribution to the column density, and hence
 *       airmass, is about 10^5 times smaller than the contribution from the
 *       lower layers, so that this assumption contributes an error of around
 *       0.001%.  If you desire more precision, you could change the upper
 *       limit on the integration in the code, but to do it properly you would
 *       have to deal with the physics of the thermosphere, which does not
 *       have quite the same properties as the lower layers (for example, the
 *       assumption "constancy of chemical composition" breaks down).
 *  CRC FORMULA FOR MU(WAVELENGTH):  UNIMPORTANT.  This means the error in
 *       calculating the index of refraction mu for a given wavelength, because
 *       the formula in the CRC Handbook is not exact.  The wavelength of
 *       observation already has almost no effect on the airmass calculation
 *       (see IMPORTANCE OF ATMOSPHERIC VARIABLES below); the error in using
 *       the CRC formula (which itself is not a bad estimator for the index of
 *       refraction, at least for reasonable temperatures and pressures) must
 *       therefore have even less of an impact.  I feel pretty confident in
 *       saying that this is a good assumption.  Note that the code does not
 *       include an additional term in the formula for mu which accounts for
 *       humidity; this is safe because of the smallness of the term and
 *       because airmass depends on mu so weakly.
 *
 *
 *  VERBOSITY OPTION
 *  ================
 *
 *  With the -v ("verbose") flag, the program also spits out the airmass
 *  estimated by the polynomial approximation formula, and the percent error
 *  in this approximation (relative to the supposed "true" airmass computed
 *  by this program).  It also computes the result and error of the simple
 *  formula X ~ sec(z), and of the airmass formula used within IRAF.  The
 *  latter formula comes from  J.A. Ball, _Algorithms for the HP-45 and HP-35_
 *  (1975 ed.) (Cambridge, MA: Center for Astrophysics), with Q assumed to be
 *  750.  (Ball's algorithm comes from C.W. Allen, _Astrophysical Quantities_
 *  (3rd ed.) (London: Athlone), p. 133 (1973).)
 *       For all but extreme zenith angles (i.e., for z < 86 degrees), the
 *  sec(z) - 1 polynomial is the best of the three approximations.  The
 *  polynomial only goes bad at extreme zenith angles because it was not
 *  designed to fit these large angles.  Above about 86 degrees, the IRAF
 *  formula is a better approximation.  (The sec(z) - 1 polynomial goes negative
 *  above z ~ 88.35, which is of course physically impossible.)  But at these
 *  extremes it is clearly best to use the airmass algorithm utilized in this
 *  program.
 *       For low to moderate zenith angles, i.e., z < 70 roughly, it seems
 *  perfectly safe to use the sec(z) - 1 polynomial to approximate the airmass
 *  unless very high accuracy is required or atmospheric conditions are highly
 *  abnormal.  The error of that approximation is about 0.09% for typical
 *  atmospheric conditions.  The sec(z) - 1 polynomial error exceeds 1% at
 *  about z = 85 degrees, again for typical conditions.  The approximation
 *  rapidly deteriorates above this z.
 *       To judge the conditions over which the sec(z) - 1 polynomial is
 *  suitable, here is a table for a few large zenith angles.  Again, this is
 *  for the default values of atmospheric variables.  Also tabulated are
 *  percent errors for the polynomial when atmospheric conditions are at the
 *  most extreme one could expect to encounter at any reasonable point on the
 *  earth's surface.  The sec(z) and IRAF formula errors are also included.
 *
 *                           PERCENT ERRORS RELATIVE TO THIS PROGRAM
 *       z (deg)    sec(z)-1  (extreme)   IRAF   (extreme)   sec(z)  (extreme)
 *       =======    ===================   ================   =================
 *        60          0.0346     0.12      0.111    0.20      0.310     0.40
 *        65          0.0537     0.18      0.169    0.29      0.474     0.60
 *        70          0.0889     0.27      0.273    0.46      0.774     0.96
 *        75          0.159      0.44      0.489    0.77      1.41      1.7
 *        80          0.280      0.90      1.04     1.7       3.16      3.8
 *        81          0.294      1.0       1.25     2.0       3.87      4.6
 *        82          0.279      1.2       1.52     2.4       4.84      5.8
 *        83          0.186      1.7       1.88     3.0       6.20      7.4
 *        84          0.106      2.5       2.37     3.9       8.20      9.8
 *        85          0.942      4.1       3.02     5.0      11.3      14.
 *        86          3.42       7.6       3.89     6.6      16.5      20.
 *        87         12.0       17.        4.99     8.9      26.2      31.
 *        88         52.0       56.        6.06    12.       47.7      56.
 *        88.1       61.8       65.        6.13    12.       51.2      60.
 *        88.2       73.9       76.        6.20    13.       55.2      64.
 *        88.3       89.2       90.        6.25    13.       59.7      70.
 *        88.4       ----       ----       6.28    13.       64.8      76.
 *        88.5       ----       ----       6.29    14.       70.6      82.
 *        88.6       ----       ----       6.28    14.       77.3      90.
 *        88.7       ----       ----       6.25    14.       85.1      99.
 *        88.8       ----       ----       6.19    14.       94.3     110
 *        88.9       ----       ----       6.09    15.      105.      120
 *        89         ----       ----       5.96    15.      118.      140
 *        89.5       ----       ----       4.58    16.      266.      300
 *        89.9       ----       ----       2.30    16.        1.46E3    1.7E3
 *        89.99      ----       ----       1.61    16.        1.50E4    1.7E4
 *        89.999     ----       ----       1.53    16.        1.50E5    1.7E5
 *
 *  As the above table shows, atmospheric variable effects become important at
 *  large zenith angles (on the order of z = 80 and higher, but perhaps z ~ 70
 *  or even lower angles if you need high accuracy or the deviation of the
 *  atmospheric variables from standard values is extreme).
 *       Incidentally, sec(z) *ALWAYS* overestimates the true airmass; that is,
 *  you could consider sec(z) an upper bound on the value of the true airmass
 *  if you needed a rough estimate.  (And, of course, the sec(z) estimate
 *  progressively gets worse with larger z.)
 *
 *
 *  IMPORTANCE OF ATMOSPHERIC VARIABLES
 *  ===================================
 *
 *  For low to moderate zenith angles, the dependence of airmass on atmospheric
 *  variables (pressure, temperature, etc.) is generally small unless you need
 *  high accuracy or your observing conditions are extreme.  For these zenith
 *  angles, you can probably safely assume the values built into the program by
 *  default: temperature = 15 C, pressure = 1013.25 millibars, altitude = mean
 *  sea level, latitude = 45 degrees, day of year = 80, wavelength = 5500
 *  Angstroms.  (Most of these aren't technically atmospheric conditions but
 *  are considered here as site-dependent variables like the local temperature
 *  and pressure.)  These values correspond in a rough way to the definition of
 *  the U.S. Standard Atmosphere.  (The correspondence is not exact because of
 *  the way we calculate the characteristics of the bottom two temperature
 *  layers.)
 *       The importance of these variables in affecting the airmass calculation
 *  is roughly: temperature most important; pressure of medium importance;
 *  seasonal effects (latitude combined with day of year) not very important;
 *  altitude and wavelength unimportant.  A note about the above sentence:
 *  it only refers to the effects of actually changing the PARAMETERS to the
 *  program; each parameter is associated with an effect, e.g., "-a" for
 *  altitude.  However, the "real-life" effect of changing certain parameters
 *  is strongly coupled to changing others.  That is, if the observer literally
 *  changes his *altitude* (for example, he moves further up a mountain), the
 *  locally observed *pressure* is going to change tremendously (after all,
 *  altimeters are in essence just barometers), as well as the observed
 *  temperature (it tends to get colder as you move up), so that the true
 *  effect of changing the observer's altitude is actually stronger than
 *  if the local temperature or pressure alone would change, despite what
 *  the first sentence says.  Similarly, if the observer moves a great distance
 *  in latitude, the local temperature will obviously change.  To make sure
 *  there isn't any confusion about this point: just feed in the actual values
 *  of the observing conditions (temperature, pressure, altitude, latitude,
 *  day of year, wavelength) and everything will be fine.
 *       The wavelength (the wavelength of light we're observing at) has such
 *  a weak effect because it only affects the airmass calculation indirectly,
 *  in the calculation of the path the light takes through the atmosphere; the
 *  majority of the wavelength-dependent effect enters by way of refraction,
 *  whose effect is to change the APPARENT zenith angle of the object, which
 *  of course by far is the main parameter (z) affecting airmass; non-refraction
 *  wavelength dependence is thus rather minor.  Altitude is minor because of
 *  the way the bottom two temperature layers are calculated.  (It only enters
 *  the algorithm at all through the calculation of the lower limit of the
 *  integration.)  Temperature is the most important observing condition
 *  because it determines the temperature profile within the two lowest layers,
 *  which is of tremendous importance in calculating the density (c.f.
 *  Equation (5).)
 *        Humidity is NOT treated in this code for determining the airmass
 *  of dry air.  This is because, while certainly the presence of water vapor
 *  (especially condensed water!) affects the extinction coefficient of air,
 *  it does not affect the airMASS very much.  That is, it does not greatly
 *  affect the molecular weight or total density, since even saturated air does
 *  not carry much water vapor by mass.  Furthermore, to calculate the dry
 *  airmass precisely, taking the presence of water vapor correctly into
 *  account, requires an accurate knowledge of humidity as a function of height
 *  in the atmosphere, and as discussed under EXTINCTION above, we don't
 *  usually have such information.  In any case humidity's affect on airmass
 *  (particularly since airmass is normalized to 1 at the zenith) has got to be
 *  very small, for any conceivable atmospheric condition besides perhaps,
 *  perhaps, the innermost winds of a hurricane, and at that point other
 *  assumptions like spherical symmetry break down anyway.  Humidity is
 *  implicitly included in the sense that refraction changes the apparent zenith
 *  angle of the object, and any good refraction algorithm includes humidity.
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*******************************************************************************
 *  Mathematical and physical constants, and constants derived from them.
 *  It's probably best to keep these as they are, since the temperature levels
 *  defined in the U.S. Standard Atmosphere were calculated using the values
 *  of the physical constants below.
 ******************************************************************************/

#define DEG2RAD (M_PI/180.)      /* factor for converting degrees to radians */
#define EQUATOR 6378178.           /* Earth's equatorial radius in meters */
#define POLAR (EQUATOR*(1.-1./298.32))  /* Earth's polar radius in meters */
#define GEOCENTGRAVCONST 3.9862216E14    /* G * mass of Earth (m^3 / s^2) */
#define STANDARDG0 9.80665               /* standard value of g0, in m/s^2 */
#define GASCONSTANT 8.314510       /* universal gas constant (R), in J/mol/K */
#define AIRKGPERMOLE 0.0289644    /* dry air mol.wt. 28.9644 * 0.001 g -> kg */
#define VAPORKGPERMOLE 0.0180153    /* water mol.wt. 18.0153 * 0.001 g -> kg */
#define AIRCONSTANT (GASCONSTANT/AIRKGPERMOLE)    /* R for dry air, J/kg/K */
#define VAPORCONSTANT (GASCONSTANT/VAPORKGPERMOLE)  /* R for vapor, J/kg/K */
#define CONSTC (AIRCONSTANT*273.15/101325.) /* constant for calculating c */
#define CONSTRHO (-(STANDARDG0)/AIRCONSTANT) /* constant for calculating rho */
#define STANDARDR02 (GEOCENTGRAVCONST/STANDARDG0)     /* GM/g0,  in m^2 */
#define CPVAPOR 1.81E3    /* typical heat capacity of water vapor, J/kg/K */
#define CWATER 4.20E3     /* typical specific heat of water, J/kg/K */
#define T0VAPOR 298.15   /* saturation vapor press. calibrated to this T (K) */
#define LATENTVAPOR 2443240. /* latent heat vap. at 25 C and 760 mm Hg (J/kg) */
#define P0VAPOR 3167.6       /* saturation vapor press. at 25 C (Pascals) */
#define VAPPRESSFACT1 ((CWATER-CPVAPOR)/VAPORCONSTANT) /* used in vaporrho() */
#define VAPPRESSFACT2 (((LATENTVAPOR/T0VAPOR)+(CWATER-CPVAPOR))/VAPORCONSTANT)
#define NTLEVELS 8  /* Number of different temperature regions in atmosphere */
#define K 5      /* constant used in polint() and qromb() */
#define IMAX 50  /* constant used in qromb() */

/*******************************************************************************
 *  Global variables (actually, constants, once they are computed in main().)
 *  The first two default values for betavec and tvec are replaced by values
 *  computed in main() (the default values come from the U.S. Standard
 *  Atmosphere, but we attempt a more "sophisticated" estimation of them based
 *  on the observed local temperature, pressure, etc.)
 ******************************************************************************/

static double cee, const6, const7, rhovec[NTLEVELS], bigrvec[(NTLEVELS+1)],
     betavec[NTLEVELS]={-0.0065, 0., 0.0010, 0.0028, 0., -0.0020, -0.0040, 0.},
     tvec[NTLEVELS]={288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 252.65,
     180.65}, relhumid=-1.;

/*******************************************************************************
 *  vaporrho - called only if "-H" was passed on the command line.
 *  Calculates the wator vapor density, given the temperature (t) and
 *  relative humidity (rh).
 ******************************************************************************/

inline double vaporrho(double t, double rh)
{
   return rh*(P0VAPOR/VAPORCONSTANT)*pow(T0VAPOR/t, VAPPRESSFACT1)*
        exp(VAPPRESSFACT2-(VAPPRESSFACT2*T0VAPOR)/t)/t;
}

inline double frho(double rdiff, double rhobottom, double tbottom, double beta)
{
   return (fabs(beta)<1.E-10) ? rhobottom*exp(CONSTRHO*rdiff/tbottom) :
        rhobottom*pow(tbottom/(tbottom+beta*rdiff), 1.-CONSTRHO/beta);
}

inline double musquare(double rho)
{
   return (3.+4.*cee*rho)/(3.-2.*cee*rho);
}

double columndensityint(double r)
{
   unsigned long ju=NTLEVELS, jm, index=0;
   double rho;

   /* figure out index, the temperature layer corresponding to the given r. */
   /* The algorithm is based on "locate" in Numerical Recipes, and handles */
   /* r values outside the "proper" range (we simply assume the bottom (top) */
   /* temp. layer continues on to negative (positive) infinity, in such cases)*/
   while (ju-index>1) {
      jm=(ju+index)>>1;
      if (r>bigrvec[jm]) index=jm;
      else ju=jm;
   }

   rho=frho(r-bigrvec[index], rhovec[index], tvec[index], betavec[index]);
   r*=r;
   return STANDARDR02*rho/(r*sqrt(1.-const6*r/musquare(rho)));
}

double vaporcolumndensityint(double r)
{
   double rho, rdiff, vaprho;

   if (r<bigrvec[1]) {   /* if we are in the troposphere */
      rdiff=r-*bigrvec;
      rho=frho(rdiff, *rhovec, *tvec, *betavec);
      vaprho=vaporrho(*tvec+*betavec*rdiff, relhumid);
   } else {   /* assume we're in the lower stratosphere */
      rdiff=r-bigrvec[1];
      rho=frho(rdiff, rhovec[1], tvec[1], betavec[1]);
      vaprho=vaporrho(tvec[1]+betavec[1]*rdiff, const7*(bigrvec[2]-r));
   }
   r*=r;
   return STANDARDR02*vaprho/(r*sqrt(1.-const6*r/musquare(rho)));
}

/*******************************************************************************
 *  trapzd - auxiliary to qromb integration routine
 *  Adapted from J. Percival, who in turn adapted it from the original
 *  Numerical Recipes routine by Press et al.
 ******************************************************************************/

inline double trapzd(double (*func)(double), double a, double b, int n)
{
   double del, sum, tnm, x;
   int j;
   static double s;
   static int it;

   if (n<=0) {
      s=0.5*(b-a)*((*func)(a)+(*func)(b));
      it=1;
   } else {
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      sum=0;
      for(j=0;j<it;j++) {
         sum+=(*func)(x);
         x+=del;
      }
      s=0.5*(s+(b-a)*sum/tnm);
      it*=2;
   }
   return s;
}

/*******************************************************************************
 *  polint - auxiliary to qromb integration routine
 *  Adapted from J. Percival, who in turn adapted it from the original
 *  Numerical Recipes routine by Press et al.
 ******************************************************************************/

inline double polint(double *xa, double *ya)
{
   double c[K], d[K], den, dif, dift, ho, hp, w, y;
   int i, m, ns;

   /* find the index of the closest table entry */
   ns=0;
   dif=fabs(xa[ns]);
   for(i=0;i<K;i++) {
      dift=fabs(xa[i]);
      if (dift<dif) {
         ns=i;
         dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
   }

   y=ya[ns--]; /* first guess */
   for (m=0;m<(K-1);m++) {
      for(i=0;i<(K-1)-m;i++) {
         ho=xa[i];
         hp=xa[i+m+1];
         w=c[i+1]-d[i];
         den=ho-hp;
         den=w/den;
         c[i]=ho*den;
         d[i]=hp*den;
      }
      y+=(2*ns+1<(K-2)-m) ? c[ns+1] : d[ns--];
   }
   return y;
}

/*******************************************************************************
 *  qromb - Romberg rule integration
 *  Adapted from J. Percival, who in turn adapted it from the original
 *  Numerical Recipes routine by Press et al.
 *  Changed the convergence criterion to fabs(ss-oldss)<eps.
 ******************************************************************************/

double qromb(double (*func)(double), double a, double b, double eps)
{
   double ss=0., oldss=0., h[(IMAX+1)], s[IMAX];
   int i;

   h[0]=1;
   for (i=0;i<IMAX;i++) {
      s[i]=trapzd(func, a, b, i);
      if (i>=(K-1)) {
         ss=polint(h+i-(K-1), s+i-(K-1));
         if (fabs((ss-oldss)/ss)<eps) return ss;
         oldss=ss;
      } else oldss=s[i];
      h[i+1]=0.25*h[i];
   }

   printf("EXCEEDED MAX ITERATIONS (%d) IN qromb\n", IMAX);
   return ss;
}


void usage(void)
{
   printf("Usage:\n\nairmass [-a <alt>] [-d <day>] [-H <h>] [-l <lat>] [-L "
        "<wavelength>] [-p <P>]\n        [-t <T>] [-v] <zenith>\n\nwhere "
        "<zenith> is the APPARENT zenith angle, in degrees, whose airmass you\n"
        "want to calculate.  Program by default uses standard atmosphere "
        "values;\nto override these use any of: \n\t"
        "-a <alt> : observer's altitude above sea level, in meters (default 0)"
        "\n\t-d <day> : day number from beginning of year (0 = Jan 1; default "
        "80)\n\t-H <h> : relative humidity at observer, in percent (default "
        "n/a)\n\t-l <lat> : observer's latitude in degrees (default +45)"
        "\n\t-L <wavelength> : the observing wavelength in Angstroms (default "
        "5500)\n\t-p <P> : local pressure, in millibars (default 1013.25)"
        "\n\t-t <T> : local temperature, in Celsius (default 15).\n"
        "Additionally, you can specify the -v option to calculate the airmass"
        "\npredicted by standard approximate formulas (sec(z)-1 polynomial, "
        "IRAF, sec(z))\nand their error relative to the airmass computed by "
        "this program.\nIf -H is used, program will also compute a *CRUDE* "
        "estimate of the water vapor\ncolumn density and airmass, using the "
        "specified <h>.\n");
   exit(1);
}

int compute_airmass(double alt, double daynum, double relhumid, double lat,
                    double lambda, double p0, double t0)
{
   char verbos=0;
   int i;
   double colz, col0, cosyearfrac, sinz, rmsl, bigrmsl, zfact, approximation,
        rjunk, coslatsq, airmass,
        z=-1., betapoly[5]={-0.0065107, -4.5403e-06,
        3.6599e-07, -2.2174e-09, 7.9392e-12}, r1poly[11]={17.204, 8.9155e-03,
        -3.6420e-03, 2.5617e-05, 2.4796e-07, -1.2774e-08, 1.3017e-10,
        2.0151e-12, -2.6985e-14, -1.0397e-16, 1.4849e-18}, hvec[(NTLEVELS+1)]=
        {0., 11000., 20000., 32000., 47000., 52000., 61000., 79000., 88743.};

   if ((z<0.) || (z>= 90.)) {
      printf("zenith angle not found or unacceptable value\n");
      usage();
   }

   sinz=sin(z*DEG2RAD);
   cee=CONSTC*(2.87566E-4+134.12/(lambda*lambda)+3.777E8*pow(lambda, -4.));
   cosyearfrac=cos((daynum-202.)*(360.*DEG2RAD/365.));
   coslatsq=cos(lat*DEG2RAD);
   coslatsq*=coslatsq;
   rmsl=sqrt(((POLAR*POLAR*POLAR*POLAR)+(EQUATOR*EQUATOR*EQUATOR*EQUATOR-
        POLAR*POLAR*POLAR*POLAR)*coslatsq)/((POLAR*POLAR)+(EQUATOR*EQUATOR-
        POLAR*POLAR)*coslatsq));
   bigrmsl=(-STANDARDR02)/rmsl;

   /* Calculate bigrvec, the bigr at the bottom of each temperature layer */
   *bigrvec=(-STANDARDR02)/(rmsl+alt);
   rjunk=r1poly[10];
   for (i=9;i>=0;i--) rjunk=rjunk*lat+((i%2) ? r1poly[i]*cosyearfrac :
        r1poly[i]);
   bigrvec[1]=(-STANDARDR02)/(rmsl+rjunk*1000.);
   for (i=2;i<(NTLEVELS+1);i++) bigrvec[i]=hvec[i]+bigrmsl;

   /* Set up our temperature profile for the troposphere/lower stratosphere */
   *betavec=betapoly[4];
   for (i=3;i>=0;i--) *betavec=*betavec*lat+((i%2) ? betapoly[i]*cosyearfrac :
        betapoly[i]);
   *betavec*=(-rmsl/bigrvec[1]);
   *tvec=t0;
   tvec[1]=t0+*betavec*(bigrvec[1]-*bigrvec);
   betavec[1]=(tvec[2]-tvec[1])/(bigrvec[2]-bigrvec[1]);

   /* Compute the density at the bottom of each temperature layer */
   *rhovec=p0/(AIRCONSTANT*t0);
   for (i=0;i<(NTLEVELS-1);i++) rhovec[i+1]=frho(bigrvec[i+1]-bigrvec[i],
        rhovec[i], tvec[i], betavec[i]);

   const6=musquare(*rhovec)*sinz*sinz/(*bigrvec*(*bigrvec));    /* for z */
   colz=qromb(columndensityint, *bigrvec, bigrvec[NTLEVELS], 1.E-8);
   const6=0.;  /* equivalent to setting z = 0. */
   col0=qromb(columndensityint, *bigrvec, bigrvec[NTLEVELS], 1.E-8);
   airmass=colz/col0;

   printf("for z=%.15g: column density=%.8g g/cm^2; zenith column density=%.8g g/cm^2; AIRMASS=%.8g\n", z,
        0.1*colz, 0.1*col0, airmass);

   /* if desired, compute error in using various airmass approximations */
   if (verbos) {
      printf("Comparisons to various approximate formulae for the airmass:\n");

      /* Do the sec(z)-1 polynomial approximation */
      zfact=1./cos(z*DEG2RAD)-1.;
      approximation=1.+zfact*(0.9981833-zfact*(0.002875+(zfact*0.0008083)));
      printf("sec(z) - 1 polynomial gives %.8g, an error of %.3g%%.\n",
           approximation, 100.*fabs(approximation-airmass)/airmass);

      /* Do the IRAF approximation (based on J.A. Ball, from Allen) */
      zfact=750.*cos(z*DEG2RAD);
      approximation=sqrt(zfact*zfact+1501.)-zfact;
      printf("IRAF formula gives %.8g, an error of %.3g%%.\n",
           approximation, 100.*fabs(approximation-airmass)/airmass);

      /* Do the simple airmass ~ sec(z) approximation */
      approximation=1./cos(z*DEG2RAD);
      printf("sec(z) gives %.8g, an error of %.3g%%.\n",
           approximation, 100.*fabs(approximation-airmass)/airmass);
   }

   /* if desired, compute the (very crude) estimate of vapor column density */
   if (relhumid>0.) {
      const7=relhumid/(bigrvec[2]-bigrvec[1]);
      const6=musquare(*rhovec)*sinz*sinz/(*bigrvec*(*bigrvec));    /* for z */
      colz=qromb(vaporcolumndensityint, *bigrvec, bigrvec[2], 1.E-8);
      const6=0.;  /* equivalent to setting z = 0. */
      col0=qromb(vaporcolumndensityint, *bigrvec, bigrvec[2], 1.E-8);
      printf("water vapor column density=%.8g g/cm^2; water vapor "
           "airmass=%.8g\n", 0.1*colz, colz/col0);
   }
   exit(0);
}
