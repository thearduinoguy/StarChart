# StarChart -- Sky chart program for visual and binocular astronomy
#
# Copyright (c) 2008 - 2010 by David A. Wallace
# Copyright (c) 2012 Walter Bender
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# CITATIONS:
#
#   The algorithms herein were mostly gleaned from:
#     Practical Astronomy with Your Calculator
#       by Peter Duffett-Smith (3rd ed.)
#
#   Keplerian Elements for Approximate Positions of the Major Planets
#     http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
#
#   Bright Star Catalog
#     http://heasarc.gsfc.nasa.gov/W3Browse/star-catalog/bsc5p.html
#
#   Dimmer Stars
#     The catalog of stars to magnitude 8 was extracted from the SAO catalog
#     of stars at
#     http://cdsarc.u-strasbg.fr/viz-bin/Cat?cat=I%2F131A&target=http
#     Which was filtered to only stars of visual magnitude 8.9 or brighter and
#     only those stars with a Henry Draper (HD / HDE) designator.
#
#   Names of Stars
#     extracted from data at http://www.alcyone.de/SIT/bsc/bsc.html and
#     http://cdsarc.u-strasbg.fr/viz-bin/Cat?IV/27 (table3.dat)
#
#   Messier Catalog
#     extracted from data at http://astro.nineplanets.org
#
#   Constellation figures -- derived from charts at
#     http://www.astro.wisc.edu/~dolan/constellations/
#     and the coordinates of the stars that the line-segments interconnect.
#
#   NGC Objects
#     The list of NGC objects came from several sources including:
#       http://www.hawastsoc.org/deepsky/bennett.html
#       http://www.wikipedia.org
#     
#
#  ACKNOWLEDGEMENTS
#
#  The author wishes to thank the following people for assistance with
#  this project:
#
#    Daniel Castilo and Oscar Mendez Laesprella for encouragement,
#    suggestions, bug reports, beta testing and Spanish translation.
#    Thanks, guys -- I couldn't have done this without you.
#
#    The owner and staff of The Java Room in Chelmsford, Massachusetts
#    for the coffee, wi-fi access and especially the live music.  Best
#    environment I've ever had for developing code in!
#
#    The members and officers of the Amateur Telescope Makers of
#    Boston who encouraged and educated me when I first discovered the
#    wonder that is our planet's night sky.
#
# -----------------------------------------------------------------------------
#
# INDEX
#
# (Line numbers given below are approximate, within 25 lines, usually.)
#
#   Line No.    Content
#   --------    --------------------------------------
#      175      Planetary Catalog data
#      275      Start of the code -- astronomical algorithms
#      475	convert az, alt to ra, dec
#      725      Definitions for the GUI controls (in toolbars, mostly)
#      775      Version and date displayed by "About" function.
#      800	Definition of pixel-to-object mapping class
#      925	Definition of object-to-pixel mapping class
#      950	Definition of Location class
#     1000      Definition for the ChartDisplay (main GUI) class
#     1075      Convert az,alt to x,y coordinates
#     1100      Convert x,y to az,alt coordinates
#     1125      Control callbacks method
#     1450	Code for "Locate" feature
#     1475      Magnification code
#     1600      Input Error reporting
#     1625      Code to draw the chart:
#               1625  The outline and field
#		1750  Plot the sky
#               1775  PlotAllStars code
#               1800  PlotAllDSOs code
#               1825  PlotAllConstellations code
#               1850  PlotAllPlanets code
#               2125  Compute Vector Of Magnified Object
#               2175  PlotMagStars code
#               2200  PlotMagDSOs code
#               2225  PlotMagPlanets code
#               2500  Draw symbol of Star
#               2525  Draw symbol of Planets (including sun and moon)
#               2675  Draw symbol of Deep Sky Objects
#     2850      Definition for the StarChartActivity class
#     3050      Read metadata
#     3075      Write metadata
#     3100      StarChart.cfg update code
#


# ================== NON-STANDARD BUT REQUIRED LOCALIZATION ===================

# The settings for the observer's coordinates
try:
  import observatory
  longitude = observatory.data[0]
  latitude = observatory.data[1]
except:
  longitude = 0.0
  latitude = 0.0

# "nonlocal_timezone_correction" allows your XO to keep one timezone
# while the program's "now" expresses time in another timezone.  This
# value is normally set to zero, but if you are traveling, you might
# find it more convenient to leave the XO's local timekeeping alone
# and simply offset the star chart's time zone so that "now" is
# correct for your current locale.
#
nonlocal_timezone_correction = 0.0

# If traveling, create "travel.py" to define alternate location and timezone.
# (If "travel.py" does not exist, we use our home observatory parameters.)

try:
  import travel
  longitude = travel.data[0]
  latitude = travel.data[1]
  nonlocal_timezone_correction = travel.data[2]
except:
  pass


# =================================== IMPORTS ================================

import pygtk
pygtk.require('2.0')
import gtk
import sys
import os
import gobject
from datetime import *
from time import sleep
from time import time as systime
from math import sin, cos, tan, asin, acos, atan, pi, sqrt, atan2
from sugar.activity import activity
from sugar.activity.activity import get_bundle_path
try:
    from sugar.graphics.toolbarbox import ToolbarBox
    _have_toolbox = True
except ImportError:
    _have_toolbox = False
if _have_toolbox:
    from sugar.activity.widgets import ActivityToolbarButton
    from sugar.graphics.toolbarbox import ToolbarButton
    from sugar.activity.widgets import StopButton
from sugar.graphics.toolbutton import ToolButton
from sugar.graphics.palette import Palette, ToolInvoker
from sugar.graphics.icon import Icon
import logging
from gettext import gettext as _

class ToggleButtonTool(ToolButton):

    def __init__(self, icon_on=None, icon_off=None, **kwargs):
      super(ToolButton, self).__init__(icon_off)

      self._icon_on=icon_on
      self._icon_off=icon_off
      self._accelerator = None
      self._tooltip = None
      self._palette_invoker = ToolInvoker()
      self._active = False

      gobject.GObject.__init__(self, **kwargs)

      self._palette_invoker.attach_tool(self)

      if icon_off:
        self.set_icon(icon_off)

      self.connect('clicked', self._click_callback)

    def _click_callback(self, button):
      self.set_active(not self._active)

    def get_active(self):
      return self._active

    def set_active(self, active):
      self._active = active
      if self._active:
        self.set_icon(self._icon_on)
      else:
        self.set_icon(self._icon_off)

# ============================= PLANETARY CATALOG DATA =======================
#
# (This is the only catalog that is coded into the program as opposed
# to being imported.)
#
# Keplarian data for planets (J2000.0 epoch coordinates, except moon
# is J1990.5)
#
# FIXME: add data for dwbar, dI and dO
#
#   planet        arg of peri   eccentricity radius (au)  inclination  RAAN
#   mean longitude  delta mean longitude
#   pname         wbar0         E            a            I0           O0
#   L               dL
planets = [
  # TRANS: http://en.wikipedia.org/wiki/Mercury_(planet)
  (_('Mercury'),    77.45779628, 0.20563593,  0.38709927,  7.00497902,  48.33076593, 252.25032350,   149472.67411175),
  # TRANS: http://en.wikipedia.org/wiki/Venus
  (_('Venus'),    131.60246718, 0.00677672,  0.72333566,  3.39467605,  76.67984255, 181.97909950,    58517.81538729),
  # TRANS: http://en.wikipedia.org/wiki/Earth
  (_('Earth'),    102.93768193, 0.01671123,  1.00000261, -0.00001531,   0.0,        100.46457166,    35999.37244981),
  # TRANS: http://en.wikipedia.org/wiki/Mars
  (_('Mars'),     -23.94362959, 0.09339410,  1.52371034,  1.84969142,  49.55953891,  -4.55343205,    19140.30268499),
  # TRANS: http://en.wikipedia.org/wiki/Jupiter
  (_('Jupiter'),   14.72847983, 0.04838624,  5.20288700,  1.30439695, 100.47390909,  34.39644051,     3034.74612775),
  # TRANS: http://en.wikipedia.org/wiki/Saturn
  (_('Saturn'),    92.59887831, 0.05386179,  9.53667594,  2.48599187, 113.66242448,  49.95424423,     1222.49362201),
  # TRANS: http://en.wikipedia.org/wiki/Uranus
  (_('Uranus'),   170.95427630, 0.04725744, 19.18916464,  0.77263783,  74.01692503, 313.23810451,      428.48202785),
  ]
# TRANS: http://en.wikipedia.org/wiki/Sun
sun = (_('Sun'),  282.93768193, 0.01671123,  1.00000261,  0.0,          0.0,        280.46457166,    35999.37244981)

#       name       mean lon      lon of peri     lon of Node
#       Inclination Eccentricity  radius (km)      Parallax   Offset of Epoch
#                  L0            P0              N0
#       I           e             a                phi0       tau
# TRANS: http://en.wikipedia.org/wiki/Moon
moon = (_('Moon'), 318.351648,   36.340410,      318.510107,      5.145396,   0.054900,     384401,          0.9507,    2447891.5)

# obliquity of J2000 epoch is 23.43928 degrees -- we need this value
# (in radians) for planetary calculations
eps = pi * 23.43928 / 180.0


# -----------------------------------------------------------------------------

# The bright star catalog is imported from stars1.py.
import stars1
star_chart = stars1.data

# Appending data from supplementary catalogs is just a matter of
# importing the catalog file and updating the primary dictionary with
# the supplement.  But I must allow for the possibility that the
# supplementary catalog does not exist, so the process needs a
# try/except block.

# A supplementary catalog of stars down to magnitude 6.8 is imported
# from stars2.py.  If this catalog is not present, also disable the
# "magnify" feature.
try:
  can_magnify = True
  import stars2
  supplement = stars2.data
except:
  can_magnify = False
  supplement = {}
star_chart.update(supplement)

# The constellations figures have their own catalog.  This catalog
# could potentially be replaced by locale-specific figures, but that
# will break much of the code relating to object locating and
# identifying, since the program wants to use the 88 IAU constellation
# names (or at least their abbreviations).
import constellations
figures = constellations.data

# Here we obtain catalog data for Deep Sky Objects.  dso1 is the catalog of the
# Messier objects.

# Again, we support supplementary catalogs.  (dso2 is planned to be a
# subset of the NGC catalog and dso3 ill be a subset of the IC
# catalog.)  For DSOs, we allow any or all of the catalogs to be
# missing.

# FIXME: Because the catalogs of DSOs duplicate each other, we will
# have the same object listed by multiple names.  We will need to
# disambiguate somehow (especially in the "Identify" feature, where it
# would be nice to be able to list ALL the names for the object we
# find).
try:
  import dso1
  dso_chart = dso1.data
except:
  dso_chart = []

try:
  import dso2
  dsupplement = dso2.data
except:
  dsupplement = []
dso_chart = dso_chart + dsupplement

try:
  import dso3
  dsupplement = dso3.data
except:
  dsupplement = []
dso_chart = dso_chart + dsupplement

# We define a couple of dictionaries Locate will need: one is used to get a
# constellation's ID from its name and the other is used to get a planet's
# declination from its name.
abbrev_from_name = {}
dec_from_planet = {}


# ============================= STATE INFORMATION =============================

# initial settings for display options
nightvision = False
invertdisplay = False
fliphorizontally = False
drawconstellations = True
limitingmagnitude = 4.0
saved_lmag = 4.0
# initial settings for time
specifytime = False
saved_specifytime = False
now = datetime.utcnow()
zoneoffset = -300


# ============================== START OF CODE ================================

# Because Python's trig is done in radians and I often need answers in
# degrees, these two functions are provided to convert between radians
# and degrees.

def dtor(a):
  return a * pi / 180.0


def rtod(r):
  return r * 180.0 / pi


# This function converts a decimal time to hours, minutes and seconds
# format and returns a time object.

def floattotime(f):
  h = int(f)
  mm = (f - h) * 60.0
  m = int(mm)
  s = (mm - m) * 60.0
  s = int(s)
  t = time(h, m, s)
  return t


# This function converts a decimal angle to degrees, minutes and
# seconds format and returns a string "dddmmss".

def floattoangle(f):
  h = int(f)
  mm = (f - h) * 60.0
  m = int(mm)
  s = (mm - m) * 60.0
  s = int(s)
  return '%(degree)03dd%(minute)02dm%(second)02ds' % \
         {'degree': h, 'minute': abs(m), 'second': abs(s)}


# Convert a degrees-minutes-seconds angle to fractional degrees.
# (This function will accept any character as separator but you must
# specify degrees, minutes and seconds, e.g.: 012d34m56s, or degrees
# and minutes, e.g.: 012d34m.)

def angletofloat(s):
  try:
    d = 0.0
    i = 0
    while ((i < 4) and (i < len(s))):
      c = s[i]
      c = c[0]
      if (c.isdigit()):
        i = i + 1
      else:
        break
    t = s[0: i]
    d = float(t)
    i = i + 1
    if (i < len(s)):
      t = s[i:]
      s = t
      i = 0
      while ((i < 3) and (i < len(s))):
        c = s[i]
        c = c[0]
        if (c.isdigit()):
          i = i + 1
        else:
          break
      t = s[0: i]
      d = d + float(t) / 60.0    
      i = i + 1
      if (i < len(s)):
        t = s[i:]
        s = t
        i = 0
        while ((i < 3) and (i < len(s))):
          c = s[i]
          c = c[0]
          if (c.isdigit()):
            i = i + 1
          else:
            break
        t = s[0: i]
        d = d + float(t) / 3600.0    
    return d
  except:
    return -1.0


# Converts a utc date-time to julian day

def jtime(dd):
  y = dd.year
  m = dd.month
  d = dd.day + dd.hour / 24.0 + dd.minute / 1440.0 + dd.second / 86400.0
  if (m < 3):
    m = m + 12
    y = y - 1
  a = int(y / 100.0)
  b = 2 - a + int(a / 4)
  c = int(365.25 * y)
  cd = int(30.6001 * (m + 1))
  t = b + c + cd + d + 1720994.5
  return t


# Convert a utc date-time to gst (Greenwitch Siderial Time)

def gst(d):
  d1 = datetime(d.year, d.month, d.day, 0, 0, 0)
  t1 = d.time()
  j = jtime(d1)
  s = j - 2451545.0
  t = s / 36525.0
  t0 = 6.697374558 + 2400.051336 * t + 0.000025862 * t * t
  while (t0 < 0.0):
    t0 = t0 + 24.0
  while (t0 > 24.0):
    t0 = t0 - 24.0
  ut = t1.hour + t1.minute / 60.0 + t1.second / 3600.0
  ut = ut * 1.002737909
  t0 = t0 + ut
  while (t0 < 0.0):
    t0 = t0 + 24.0
  while (t0 > 24.0):
    t0 = t0 - 24.0
  return t0


# Convert gst to lst (local siderial time)

def lst(g):
  l = longitude / 15.0
  lst = g + l
  while (lst < 0.0):
    lst = lst + 24.0
  while (lst > 24.0):
    lst = lst - 24.0
  return lst


# convert RA to hour angle

def hrangl(ra, ut):
  g = gst(ut)
  l = lst(g)
  h = l - ra
  while (h < 0.0):
    h = h + 24.0
  while (h > 24.0):
    h = h - 24.0
  return h


# Adjust J2000 equatorial coordinates for precession.

def epochpolartonow(polar, ut):
  dec = polar[1]
  ra = polar[0]
# this method is an approximation which fails when the object is near
# the celestial pole.  So we return the epoch's ra -- a better answer
# than what the algorithm would produce
  if (dec > 88.0):
     return (ra, dec)
  if (dec < -88.0):
     return (ra, dec)
# dec is in tolerance.
  n = (ut.year - 2000.0) + (ut.month - 1.0) / \
      12.0 + 15.0 / 360.0
  ra = ra * 15.0
  ra1 = ra + ((3.0742 + \
        1.33589 * sin(dtor(ra)) * \
        tan(dtor(dec))) / 3600.0)  * n
  dec1 = dec + (20.0383 * cos(dtor(ra))) / 3600.0 * n
  return (ra1 / 15.0, dec1)


# Convert equatorial coordinates to azimuth, altitude

def polartoazalt(polar, ut):
  ra = polar[0]
  dec = dtor(polar[1])
  h = dtor(hrangl(ra, ut) * 15.0)
  l = dtor(latitude)
  a = asin(sin(dec) * sin(l) + cos(dec) * cos(l) * cos(h))
  az = acos((sin(dec) - sin(l) * sin(a)) / (cos(l) * cos(a)))
  sh = sin(h)
  if (sh > 0.0):
    az = 2.0 * pi - az
  return (rtod(az), rtod(a))


# Convert azimuth, altitude to equatorial coordinates

def azalttopolar(azalt, ut):
  az = dtor(azalt[0])
  a = dtor(azalt[1])
  h = dtor(lst(gst(ut)) * 15.0)
  l = dtor(latitude)
  dec = asin(sin(a) * sin(l) + cos(a) * cos(l) * cos(az))
  hp = acos((sin(a) - sin(l) * sin(dec)) / (cos(l) * cos(dec)))
  if (sin(az) > 0):
    hp = 2.0 * pi - hp
  ra = h - hp
  if (ra < 0.0):
    ra = ra + 2.0 * pi
  return (rtod(ra), rtod(dec))


# solve M = E - (57.29578 * e) * sin(E) for E
# (not currently used -- the planetary results are accurate enough without it.)

def eccentric_anomaly(M, e):
  dE = .05
  estar = 180.0 * e / pi
# Iterate dM = M - En - estar * sin(En); dE = dM / (1 - e * cos(En));
# En = En + dE until abs(dE) <= 0.001 note that M and estar are in
# degrees, as is E and dE so make the appropriate conversion to
# radians before doing trig.
  En = M
  while (abs(dE) > .001):
    dM = M - (En - estar * sin(dtor(En)))
    dE = dM / (1 - e * cos(dtor(En)))
    En = En + dE
  return En

  
# Determine the number of minutes from GMT specified by the zone offset string.

def parse_zone_offset(s):
  try:
    bneg = False
    oh = 0
    om = 0
    i = 0
    if (s[0] == '-'):
      bneg = True
      i = 1
    elif (s[0] == '+'):
      i = 1
    j = i
    while (s[j].isdigit()):
      j = j + 1
    if (j == i):
      return (0, 0)
    oh = int(s[i: j])
    i = j
    if (s[i] != ':'):
      return (0, 0)
    i = i + 1
    om = int(s[i:])
    if (bneg):
      oh = -oh
      om = -om
    return (oh, om)
  except:
    return (0, 0)


# Determine the datetime value represented by the timestamp string.

def parse_timestamp(s):
  try:
    if ((s[4] != '/') or (s[7] != '/') or (s[10] != ',') or (s[13] != ':')):
      return (2000, 1, 1, 0, 0)
    Y = int(s[0: 4])
    M = int(s[5: 7])
    D = int(s[8: 10])
    h = int(s[11: 13])
    m = int(s[14:])
    return (Y, M, D, h, m)
  except:
    return (2000, 1, 1, 0, 0)


# Get the GMT datetime value from the timestamp and offset strings

def get_time_and_UTC_offset(timestr, offsetstr):
  class TZ(tzinfo):
    def __init__(self, offset, name):
        self.__offset = timedelta(minutes = offset)
        self.__name = name

    def utcoffset(self, dt):
        return self.__offset

    def tzname(self, dt):
        return self.__name

    def dst(self, dt):
        return timedelta(0)

# Parse the timestamp string and the UTC offset string into year, month, day,
# hour and minute and hour and minute, respectively.

  (sy, sM, sd, sh, sm) = parse_timestamp(timestr)
  (oh, om) = parse_zone_offset(offsetstr)

# now convert the zone offset to a timezone object
             
  tzo = oh * 60 + om
  tz = TZ(tzo, '')

# using the parsed time and timezone object, construct a datetime object that
# represents the local time.
             
  lt = datetime(sy, sM, sd, sh, sm, 0, 0, tz)

# finally, subtract the zone object from the local time to get GMT.

  gt = lt - tz.utcoffset(0)
  return gt


# Using the current local time, the nonlocal_timezone_correction (if any) and
# the GMT time, set the text in entry3 and entry4 to the "real" local time and
# the "real" zone offset, respectively

def set_time_and_UTC_offset():
# convert nonlocal_timezone_correction to a time displacement
  to = timedelta(0, 60 * nonlocal_timezone_correction)
  gt = datetime.utcnow()
  lt = datetime.fromtimestamp(systime())
  dt = lt - gt + to
  tt = lt
  tt = tt + to
  lts = tt.strftime('%Y/%m/%d,%H:%M')
  dth = dt.days * 24 + int(dt.seconds / 3600)
  dtm = (dt.seconds / 60) % 60
  utos = '%d:%02d' % (dth, dtm)
  return (lts, utos)


def syntax_check_time():
  try:
    s = entry3.get_text()
# this string must conform to YYYY/MM/DD,HH:MM format and be a valid
# date and time.
    if (len(s) != 16):
      return False
    if ((s[4] != '/') or (s[7] != '/') or (s[10] != ',') or (s[13] != ':')):
      return False
    z = s[0: 4]
    if (not z.isdigit()):
      return False
    y = int(z)
    z = s[5: 7]
    if (not z.isdigit()):
      return False
    x = int(z)
    if ((x < 1) or (x > 12)):
      return False
    z = s[8: 10]
    if (not z.isdigit()):
      return False
    d = int(z)
    if ((d < 1) or (d > 31)):
      return False
    if ((x == 1) or (x == 3) or (x == 5) or (x == 7) or (x == 8) or (x == 10) or (x == 12)):
      pass
    elif (x == 2):
      if (d > 29):
        return False
      elif (((y % 4) != 0) and (d > 28)):
        return False
      elif (((y % 100) == 0) and (d > 28)):
        if ((y % 400) != 0):
          return False
    elif (d > 30):
      return False
    z = s[11: 13]
    if (not z.isdigit()):
      return False
    h = int(z)
    if (h > 23):
      return False
    z = s[14]
    if (not z.isdigit()):
      return False
    m = int(z)
    if (m > 59):
      return False
    return True
  except:
    return False


def syntax_check_zone():
  try:
    s = entry4.get_text()
# this string must conform to [+|-]HH:MM and must be a valid time between -14:00 and 14:00.
    l = len(s)
    if (l < 4):
      return False
    if ((s[0] == '+') or (s[0] == '-')):
      p = 1
    else:
      p = 0
    if ((l - p) < 4):
      return False
    if (not s[p].isdigit()):
      return False
    if (s[p+1].isdigit()):
      n = 1
    else:
      n = 0
    if ((l - p - n) < 4):
      return False
    if (s[1 + p + n] != ':'):
      return False
    if ((not s[p + 2 + n].isdigit()) or (not s[p + 3 + n].isdigit())):
      return False
    z = s[p: p + 1 + n]
    if (int(z) > 14):
      return False
    z = s[p + 2 + n : p + n + 4]
    if (int(z) > 59):
      return False
    return True
  except:
    return False


def get_planet_index(name):
# return the index number corresponding to the name of the planet asked for.
# (Mercury = 0, Venus = 1, Moon = 2, Mars = 3, Jupiter = 4, Saturn = 5,
# Uranus = 6, Sun = 7)
  (pname, L0, P0, N0, I, e, a, phi0, tau) = moon
  if (name == pname):
    return 2
  (pname, wbar, e, a, I, O, L0, dL) = sun
  if (name == pname):
    return 7
  for i in range(len(planets)):
    (pname, wbar, e, a, I, O, L0, dL) = planets[i]
    if (name == pname):
      return i
  return 6


# -----------------------------------------------------------------------------
#
# These controls affect the program's state-variables and must be global or we
# can't set or retrieve their values everywhere necessary:
#
# controls on menubar1 (_('what')):
fullscreen = ToolButton('view-fullscreen')
button1 = ToggleButtonTool(icon_off='night-off', icon_on='night-on')
button1.set_tooltip(_('Night Vision'))
button2 = ToggleButtonTool(icon_off='invert-off', icon_on='invert-on')
button2.set_tooltip(_('Invert Display'))
button3 = ToggleButtonTool(icon_off='left-right', icon_on='right-left')
button3.set_tooltip(_('Flip L/R'))
button4 = ToggleButtonTool(icon_off='constellations-off',
                           icon_on='constellations-on')
button4.set_tooltip(_('Draw Constellations'))
container2 = gtk.Table(columns=6, rows=1)
# TRANS: http://en.wikipedia.org/wiki/Magnitude_(astronomy)
label6 = gtk.Label(_('Mag:'))
rb7 = gtk.RadioButton(None, _('1'))
rb8 = gtk.RadioButton(rb7, _('2'))
rb9 = gtk.RadioButton(rb7, _('3'))
rb10 = gtk.RadioButton(rb7, _('4'))
rb11 = gtk.RadioButton(rb7, _('5'))
rb12 = gtk.RadioButton(rb7, _('6'))
# controls on menubar2 (_('where')):
container3 = gtk.VBox()
container4 = gtk.VBox()
# TRANS: http://en.wikipedia.org/wiki/Longitude
label1 = gtk.Label(_('Longitude:'))
entry1 = gtk.Entry()
entry1.set_width_chars(10)
# TRANS: http://en.wikipedia.org/wiki/East
rb1 = gtk.RadioButton(None, _('E'))
# TRANS: http://en.wikipedia.org/wiki/West
rb2 = gtk.RadioButton(rb1, _('W'))
# TRANS: http://en.wikipedia.org/wiki/Latitude
label2 = gtk.Label(_('Latitude:'))
entry2 = gtk.Entry()
entry2.set_width_chars(10)
# TRANS: http://en.wikipedia.org/wiki/North
rb3 = gtk.RadioButton(None, _('N'))
# TRANS: http://en.wikipedia.org/wiki/South
rb4 = gtk.RadioButton(rb3, _('S'))
icon = Icon(icon_name='dialog-ok')
button5 = gtk.Button()
button5.set_image(icon)
icon.show()
button5.set_label(_('Ok'))
button5.show()
button51 = ToolButton('home')
button51.set_tooltip(_('Make home'))
button51.show()
# controls on menubar3 (_('when')):
rb5 = gtk.RadioButton(None, _('Now'))
rb6 = gtk.RadioButton(rb5, _('Specify'))
label4 = gtk.Label(_('Time:'))
entry3 = gtk.Entry()
entry3.set_width_chars(16)
label5 = gtk.Label(_('Offset:'))
entry4 = gtk.Entry()
entry4.set_width_chars(7)
icon = Icon(icon_name='dialog-ok')
button6 = gtk.Button()
button6.set_image(icon)
icon.show()
button6.set_label(_('Ok'))
button6.show()
# controls on menubar4 (_('Locate')):
labell1 = gtk.Label(_('Object type:'))
objtypecb = gtk.combo_box_new_text()
planetscb = gtk.combo_box_new_text()
constscb = gtk.combo_box_new_text()
starscb = gtk.combo_box_new_text()
container0 = gtk.HBox()
container1 = gtk.VBox()
dsoscb = gtk.combo_box_new_text()
# controls on last menubar (_('About')):
# labela1 = gtk.Label(_('Version 2.0 (build 115) of 2010.04.21.1530 UT'))
# labela2 = gtk.Label(' ')
labela3 = gtk.Label(_('See http://wiki.laptop.org/go/StarChart for help.'))
labela4 = gtk.Label(' ')

# -------------------------------------------------------------------------------

# Set control states and values from the state variables.

def initialize_controls():
  button1.set_active(nightvision)
  button2.set_active(invertdisplay)
  button3.set_active(fliphorizontally)
  button4.set_active(drawconstellations)
  rb12.set_active(limitingmagnitude >= 6.0)
  rb11.set_active((limitingmagnitude >= 5.0) and (limitingmagnitude < 6.0))
  rb10.set_active((limitingmagnitude >= 4.0) and (limitingmagnitude < 5.0))
  rb9.set_active((limitingmagnitude >= 3.0) and (limitingmagnitude < 4.0))
  rb8.set_active((limitingmagnitude >= 2.0) and (limitingmagnitude < 3.0))
  rb7.set_active((limitingmagnitude >= 1.0) and (limitingmagnitude < 2.0))
#  rb7.set_active(limitingmagnitude < 1.0)
  entry2.set_text(floattoangle(abs(latitude)))
  rb4.set_active(latitude < 0.0)
  rb3.set_active(latitude >= 0.0)
  entry1.set_text(floattoangle(abs(longitude)))
  rb2.set_active(longitude < 0.0)
  rb1.set_active(longitude >= 0.0)
  rb6.set_active(specifytime)


# ========================= PixelsToObjectMap Class ============================
#
# This code is central to implementing the  "Identify" feature.

class PixelsToObjectMap():
  def __init__(self, context):
    self.data = {}
    self.context = context


  def clear(self):
    self.data.clear()


  def add(self, x, y, type, name):
    self.data[(x, y)] = (type, name)


  def identify(self, x, y):
# find the brightest object nearest to (x,y) in the map. report its details.
# (search to a distance of 16 pixels before giving up.)  Will favor brightest
# over nearest, because of the way the search algorithm is coded.  (which may be a bug)
    if self.found(x, y):
# The user got lucky and hit the target dead on.
      (type, name) = self.data[(x, y)]
      self.report(type, name)
    else:
      if (self.context.chart.magnifying):
        w = 24
      else:
        w = 8
# Search an area (2*w) pixels square for an object.
      for i in range(1, w + 1):
        xm = x - i
        xx = x + i
        ym = y - i
        yx = y + i
        (xf, yf) = self.search_square(xm, xx, ym, yx)
        if (xf > xx):
          type = ''
          name = ''
        else:
          (type, name) = self.data[(xf, yf)]
        if (type != ''):
          self.report(type, name)
        else:
          self.report(_('There is no object here.'), '')


  def get_count(self):
    return len(self.data.keys())


  def found(self, x, y):
    return (x, y) in self.data


  def search_square(self, left, right, top, bottom):
# Find all objects on the search square and return the brightest in
# event of their being more than one.  Return the object's (x,y).
#
# An extended object is always be favored over a point object if the
# square being searched includes the coordinates of the extended
# object, since the extended object is drawn over the point objects
# and therefore hides them.
    lfx = right + 1
    lfy = bottom + 1
    lfmag = 100
# For the Sun, Moon and planets we will define their magnitude to be a
# rough average: by this definition, the Sun is -27, the Moon is -6.0,
# Venus is -3.5, Mercury is -1.5, Jupiter is -0.5, Saturn is +0.5,
# Mars is +1.0 and Uranus is +5.5.  These values are a bit arbitrary,
# especially as planets can vary quite a lot, but are defined this way
# in order to ensure we have an unambiguous brightness relationship,
# irrespective of the planet's position in its orbit or phase of
# illumination.
    pmag = [ -1.5, -3.5, -6.0, 1.0, -0.5, 0.5, 5.5, -27.0 ]
    for x in range(left, right + 1):
      for y in range(top, bottom + 1):
        if self.found(x, y):
          (type, name) = self.data[(x, y)]
          if (type != 'star') and (type != 'planet'):
# Extended objects are special-case.
            if (type == 'dso-messier') or (type == 'dso-ngc'):
# Simply return the coordinates at which the object was found -- we
# don't actually need the object's catalog data.
              return (x, y)
# This object is not an extended object and therefore we will care
# about the object's brightness.  So see if this object is brighter
# than the last one found, if any.  If so, replace the last-found
# object with this one; if not, ignore this object.
          elif (type == 'planet'):
            mag = pmag[get_planet_index(name)]
          elif (type == 'star'):
# Locate the catalog entry for the star and obtain its brightness.
            (ra, dec, mag, cid) = star_chart[name]
          else:
            mag = 99 # placeholder -- object type is unknown.
          if (mag < lfmag):
            lfx = x
            lfy = y
            lfmag = mag
    return (lfx, lfy)


  def report(self, type, name):
# FIXME: for now, we only report the object's name and type.
# Eventually, we will report more and that will mean a more elaborate,
# multi-line "identification" object with multiple fields.
    if (name == ''):
      self.context.identifyobject.set_label(type)
    else:
      if (type == 'dso-messier') or (type == 'dso-ngc'):
        type = _('deep-sky object')
      elif (type == 'planet'):
        # TRANS: http://en.wikipedia.org/wiki/Planet
        type = _('planet')
        if (name == _('Sun')):
          self.context.identifyobject.set_label(_('Object is: ') + name)
          return
        if (name == _('Moon')):
          self.context.identifyobject.set_label(_('Object is: ') + name)
          return
      elif (type == 'star'):
        type = _('star')
      self.context.identifyobject.set_label(_('Object is ') + type + ': ' + name)



# ========================= ObjectToPixelsMap Class ============================
#
# This code is central to implementing the "Locate" feature.

class ObjectToPixelsMap():
  def __init__(self, context):
    self.data = {}
    self.context = context


  def clear(self):
    self.data.clear()


  def add(self, type, name, x, y):
    self.data[(type, name)] = (x, y)

  def locate(self, type, name):
# Return the coordinates for the specified object type and name.
# Returns (-1, -1) if the object is not in the map (i.e.: was not plotted).
    if self.found(type, name):
      return self.data[(type, name)]
    else:
      return (-1, -1)


  def found(self, type, name):
    return (type, name) in self.data
    

# ============================= Location Class ================================
#
# This class holds the x,y pixel coordinates of the located object so that the
# plotchart() method of the ChartDisplay class can plot a cross on that object.

class Location():
  def __init__(self, context):
    self.data = [False, 0, 0]
    self.context = context


  def clear(self):
    self.data = [False, 0, 0]


  def set(self, x, y):
    self.data = [True, x, y]


  def is_set(self):
    return self.data[0]

  def plot_cross(self):
    if (self.is_set()):

#  Draw a cross using heavy green lines, centered on the object with
#  arms which are 50 pixels long.  FIXME: As the display updates, the
#  cross is always at x,y but the highlighted object eventually drifts
#  westward.  The code needs to compensate for that somehow.  (Or else
#  this is a feature, since it shows how the sky moves over time.)

      x = self.data[1] + 1
      y = self.data[2] + 1
      self.context.gc.set_foreground(self.context.colors[4])
      self.context.gc.set_line_attributes(5, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                gtk.gdk.JOIN_MITER)
      self.context.window.draw_line(self.context.gc, x, y - 25, x, y + 25)
      self.context.window.draw_line(self.context.gc, x - 25, y, x + 25, y)
      self.context.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                gtk.gdk.JOIN_MITER)
      self.context.gc.set_foreground(self.context.colors[1])
    else:
      pass


# ============================== ChartDisplay Object ============================

class ChartDisplay(gtk.DrawingArea):
  def __init__(self, context):
    super(ChartDisplay, self).__init__()
    self.context = context
    self.colors = {}
    self.canplot = False
    self.pangolayout = self.create_pango_layout('')
    self.add_events(gtk.gdk.BUTTON_PRESS_MASK | gtk.gdk.BUTTON1_MOTION_MASK |
	gtk.gdk.BUTTON2_MOTION_MASK)
    self.connect('button_press_event', self.pressing)
    self.magnifying = False
    self.mag_center = [0, 0]
    if (not specifytime):
      gobject.timeout_add(60000, self.timer1_cb)

# The pmap object maintains a map between the chart display's pixel coordinates
# (x, y) and a visible object (type, name) for the "Identify" feature.
# The omap object maintains a map between the visible object (type, name) and
# the chart display's pixel coordinates (x, y) for the "Locate" feature.

    self.pmap = PixelsToObjectMap(context)
    self.omap = ObjectToPixelsMap(context)

# The Location object stores the result of a Locate operation.

    self.location = Location(self)

  def area_expose_cb(self, area, event):

# Determine the area we can draw upon and adjust the chart accordingly.

    rect = self.get_allocation()
    self.screensize = (rect.width, rect.height)
    self.margin = 40
    self.diameter = min(self.screensize[0], self.screensize[1]) - \
                    2 * self.margin
    self.xoffset = (self.screensize[0] - self.diameter) / 2 - self.margin
    self.yoffset = (self.screensize[1] - self.diameter) / 2 - self.margin

# Establish color selections (need only do this once).

    if (len(self.colors) == 0):
      self.gc = self.style.fg_gc[gtk.STATE_NORMAL]
      self.colormap = self.gc.get_colormap()
      self.colors[0] = self.colormap.alloc_color('white')
      self.colors[1] = self.colormap.alloc_color('black')
      self.colors[2] = self.colormap.alloc_color('red')
      self.colors[3] = self.colormap.alloc_color('gray')
      self.colors[4] = self.colormap.alloc_color('green')
      self.canplot = True
    self.plotchart()


# Catch the one-minute timer interrupt

  def timer1_cb(self):
# do not redraw the chart if we're not advancing time.
    if (not specifytime):
      self.plotchart()
    return True


# Convert az, alt (polar) coordinates to pixel position (x, y)
 
  def azalttoxy(self, polar):
    az = dtor(polar[0])
    alt = polar[1]

# radius is proportional to 90-alt

    r = self.diameter / 2.0 - alt * (self.diameter / 2.0 / 90.0)

# convert r and az to x and y
# draw the chart so east is on the right by default.
# When flipped, east is on the left, like a map.

    if (fliphorizontally):
      x = int(self.diameter / 2.0 + r * sin(az))
    else:
      x = int(self.diameter / 2.0 - r * sin(az))
    y = int(self.diameter / 2.0 - r * cos(az))

# Return the pixel coordinate relative to the center of the plot.

    return (x,y)


# Convert pixel position (x, y) to horizon (az, alt) coordinates

  def xytohorizon(self, (x, y)):
    dy = (self.margin - 2 + self.yoffset + self.diameter / 2.0) - y
    dx = (self.margin - 2 + self.xoffset + self.diameter / 2.0) - x
# Compensate for which way the chart is flipped.
    if (fliphorizontally):
      dx = -dx
    r = sqrt(dx * dx + dy * dy)
# prevent divide-by-zero:
    if (r == 0):
      return (0, 90.0)
    az = asin(dx / r)
# determine proper quadrant for az -- origin is north and goes clockwise
    if (dx >= 0) and (dy >= 0): # north to east
      pass
    elif (dx >= 0) and (dy < 0): # east to south
      az = pi - az
    elif (dx < 0) and (dy < 0): # south to west
      az = pi - az
    elif (dx < 0) and (dy >= 0): # west to north
      az = 2.0 * pi + az
# alt is inversely proportional to radius
    alt = 90.0 - 90.0 * r / (self.diameter / 2.0)
    return (rtod(az), alt)

  def callback(self, widget, data=None):

# Handle control changes here.

    global nightvision
    global invertdisplay
    global fliphorizontally
    global drawconstellations
    global limitingmagnitude
    global saved_lmag
    global longitude
    global latitude
    global specifytime

    if (data == None):
      return False
    elif (data == 'night vision'):
      nightvision = button1.get_active()
      self.plotchart()
      return False
    elif (data == 'invert display'):
      invertdisplay = button2.get_active()
      self.plotchart()
      return False
    elif (data == 'flip horizontally'):
      fliphorizontally = button3.get_active()
      self.plotchart()
      return False
    elif (data == 'draw constellations'):
      drawconstellations = button4.get_active()
      self.plotchart()
      return False
    elif (data == 'home location set'):
      s = entry1.get_text()
      if (s.find('m') >= 0):
        lon = angletofloat(s)
      else:
        try:
          lon = float(s)
        except:
          lon = -1.0
      if ((lon < 0.0) or (lon > 180.0)):
        self.lon_error()
        entry1.set_text(floattoangle(abs(longitude)))
        if (longitude < 0.0):
          rb2.clicked()
        else:
          rb1.clicked()
        return True
      else:
        longitude = lon
        entry1.set_text(floattoangle(abs(longitude)))
        if (rb2.get_active()):
          longitude = -longitude
      s = entry2.get_text()
      if (s.find('m') >= 0):
        lat = angletofloat(s)
      else:
        try:
          lat = float(s)
        except:
          lat = -1.0
      if ((lat < 0.0) or (lat > 90.0)):
        self.lat_error()
        entry2.set_text(floattoangle(abs(latitude)))
        if (latitude < 0.0):
          rb4.clicked()
        else:
          rb3.clicked()
        return True
      else:
        latitude = lat
        entry2.set_text(floattoangle(abs(latitude)))
        if (rb4.get_active()):
          latitude = -latitude
      self.location.clear()
      dsoscb.set_active(-1)
      starscb.set_active(-1)
      planetscb.set_active(-1)
      self.context.identifyobject.set_label('')
      self.context.update_config()
      self.plotchart()
      return False
    elif (data == 'location change'):
      s = entry1.get_text()
      if (s.find('m') >= 0):
        lon = angletofloat(s)
      else:
        try:
          lon = float(s)
        except:
          lon = -1.0
      if ((lon < 0.0) or (lon > 180.0)):
        self.lon_error()
        entry1.set_text(floattoangle(abs(longitude)))
        if (longitude < 0.0):
          rb2.clicked()
        else:
          rb1.clicked()
        return True
      else:
        longitude = lon
        entry1.set_text(floattoangle(abs(longitude)))
        if (rb2.get_active()):
          longitude = -longitude
      s = entry2.get_text()
      if (s.find('m') >= 0):
        lat = angletofloat(s)
      else:
        try:
          lat = float(s)
        except:
          lat = -1.0
      if ((lat < 0.0) or (lat > 90.0)):
        self.lat_error()
        entry2.set_text(floattoangle(abs(latitude)))
        if (latitude < 0.0):
          rb4.clicked()
        else:
          rb3.clicked()
        return True
      else:
        latitude = lat
        entry2.set_text(floattoangle(abs(latitude)))
        if (rb4.get_active()):
          latitude = -latitude
      self.location.clear()
      dsoscb.set_active(-1)
      starscb.set_active(-1)
      planetscb.set_active(-1)
      self.context.identifyobject.set_label('')
      self.plotchart()
      return False
    elif (data == 'user time'):
      specifytime = True
      saved_specifytime = True
      return True
    elif (data == 'now time'):
      specifytime = False
      saved_specifytime = False
      now = datetime.utcnow()
      (tstr, ostr) = set_time_and_UTC_offset()
      entry3.set_text(tstr)
      entry4.set_text(ostr)
      self.plotchart()
      return True
    elif (data == 'time change'):
      specifytime = rb6.get_active()
      saved_specifytime = specifytime
      if (specifytime):
        if (syntax_check_time() == False):
          self.time_error()
          specifytime = False
          saved_specifytime = False
          now = datetime.utcnow()
          (tstr, ostr) = set_time_and_UTC_offset()
          entry3.set_text(tstr)
          entry4.set_text(ostr)
          return True
        else:
          if (syntax_check_zone() == False):
            self.zone_error()
            specifytime = False
            saved_specifytime = False
            now = datetime.utcnow()
            (tstr, ostr) = set_time_and_UTC_offset()
            entry3.set_text(tstr)
            entry4.set_text(ostr)
            return True
      self.location.clear()
      dsoscb.set_active(-1)
      starscb.set_active(-1)
      planetscb.set_active(-1)
      self.context.identifyobject.set_label('')
      self.plotchart()
      return False
    elif (data == 'rb7 clicked'):
      limitingmagnitude = 1.0
      saved_lmag = 1.0
      self.plotchart()
      return False
    elif (data == 'rb8 clicked'):
      limitingmagnitude = 2.0
      saved_lmag = 2.0
      self.plotchart()
      return False
    elif (data == 'rb9 clicked'):
      limitingmagnitude = 3.0
      saved_lmag = 3.0
      self.plotchart()
      return False
    elif (data == 'rb10 clicked'):
      limitingmagnitude = 4.0
      saved_lmag = 4.0
      self.plotchart()
      return False
    elif (data == 'rb11 clicked'):
      limitingmagnitude = 5.0
      saved_lmag = 5.0
      self.plotchart()
      return False
    elif (data == 'rb12 clicked'):
      limitingmagnitude = 6.0
      saved_lmag = 6.0
      self.plotchart()
      return False
    elif (data == 'objtype sel'):
# Get selection and expose object selector control(s).
      selstr = objtypecb.get_active_text()
      if (selstr == _('Planets')):
        for i in reversed(range(len(container0.get_children()))):
          container0.remove(container0.get_children()[i])
        container0.add(planetscb)
        planetscb.show()
        container0.show()
        self.location.clear()
        dsoscb.set_active(-1)
        starscb.set_active(-1)
        planetscb.set_active(-1)
        self.context.identifyobject.set_label('')
      elif (selstr == _('Stars by Constellation')):
        for i in reversed(range(len(container0.get_children()))):
          container0.remove(container0.get_children()[i])
        container0.add(constscb)
        container0.add(starscb)
        starscb.get_model().clear()
        starscb.show()
        constscb.show()
        container0.show()
        self.location.clear()
        dsoscb.set_active(-1)
        starscb.set_active(-1)
        planetscb.set_active(-1)
        self.context.identifyobject.set_label('')
      elif (selstr == _('Brightest Stars')):
        for i in reversed(range(len(container0.get_children()))):
          container0.remove(container0.get_children()[i])
        container0.add(starscb)
        starscb.get_model().clear()
# Load the combobox with the names of stars whose magnitude is +1.50
# or brighter
        names = []
        for name, (ra, dec, mag, cid) in star_chart.iteritems():
          if (mag <= 1.50):
            names = names + [name]
        for name in sorted(names):
          starscb.append_text(name)
        starscb.show()
        constscb.show()
        container0.show()
        self.location.clear()
        dsoscb.set_active(-1)
        starscb.set_active(-1)
        planetscb.set_active(-1)
        self.context.identifyobject.set_label('')
      elif (selstr == _('Deep-sky Objects')):
        for i in reversed(range(len(container0.get_children()))):
          container0.remove(container0.get_children()[i])
        container0.add(dsoscb)
        dsoscb.show()
        container0.show()
        self.location.clear()
        dsoscb.set_active(-1)
        starscb.set_active(-1)
        planetscb.set_active(-1)
        self.context.identifyobject.set_label('')
      else:
        self.plotchart()
        self.location.clear()
        dsoscb.set_active(-1)
        starscb.set_active(-1)
        planetscb.set_active(-1)
        self.context.identifyobject.set_label('')
      return False
    elif (data == 'constellation sel'):
      csel = constscb.get_active_text()
      if (csel == None):
        pass
      else:
        const_id = abbrev_from_name[csel]
        starscb.get_model().clear()
# Load the stars combobox with the names of all stars having this
# constellation ID.
        names = []
        for name, (ra, dec, mag, cid) in star_chart.iteritems():
          if (cid == const_id):
            names = names + [name]
        for name in sorted(names):
          starscb.append_text(name)
        starscb.show()
      return False
    elif (data == 'planet sel'):
      selstr = planetscb.get_active_text()
      self.locate('planet', selstr)
      return False
    elif (data == 'star sel'):
      selstr = starscb.get_active_text()
      self.locate('star', selstr)
      return False
    elif (data == 'dso sel'):
      selstr = dsoscb.get_active_text()
      self.locate('dso', selstr)
      return False


  def locate(self, type, name):
    if (name == None):
      return
    if (type == 'planet'):
# Selected a planet.  Get its DEC.  If visible, highlight it;
# otherwise, call not_located.
      dec = dec_from_planet[name]
    elif (type == 'star'):
# Selected a star.  Determine which and get its DEC.  If visible,
# highlight it; otherwise, call not_located.
      (ra, dec, mag, cid) = star_chart[name]
    elif (type == 'dso'):
# Selected a DSO.  Determine which and get its DEC.  If visible,
# highlight it; otherwise, call not_located.
      for i in range(len(dso_chart)):
        (nM, strCon, ra, dec, mag, majA, minA, posA, strT, strN) = dso_chart[i]
        if (strN == ''):
          if (name == nM):
            break
        elif (name == strN + ' (' + nM + ')'):
          break
      if (nM[0:1] == 'M'):
        type = 'dso-messier'
      elif (nM[0:1] == 'N'):
        type = 'dso-ngc'
      else:
        type = 'dso'
    else:
# Invalid object type selection.  Ignore.
      return
    (x, y) = self.omap.locate(type, name)
    if ((x < 0) or (y < 0)):
      self.not_located(self.never_visible(dec))
    else:
      self.context.identifyobject.set_label('')
      self.located(x, y)


  def magnify(self, x, y):
    global limitingmagnitude
    global specifytime
    global saved_specifytime
    if (can_magnify):
# Toggle magnification state.
      if (self.magnifying):
        self.magnifying = False
        limitingmagnitude = saved_lmag
        specifytime = saved_specifytime
      else:
        self.magnifying = True
        limitingmagnitude = 9.0
        saved_specifytime = specifytime
        specifytime = True
# Save x,y.  Re-plot.
        self.mag_center[0] = x
        self.mag_center[1] = y
#      rb6.set_active(specifytime)
      self.adjust_toolbars()
      self.context.identifyobject.set_label('')
      self.location.clear()
      dsoscb.set_active(-1)
      starscb.set_active(-1)
      planetscb.set_active(-1)
      self.plotchart()
    else:
      pass


  def adjust_toolbars(self):
# Hide the "When" and "Where" toolbars when magnifying and show them
# when not.  Hide the magnitude radio buttons and the draw
# constellations button when magnifying and show these controls when
# not.

    if (self.magnifying):
      button4.hide()
      label6.hide()
      rb7.hide()
      rb8.hide()
      rb9.hide()
      rb10.hide()
      rb11.hide()
      rb12.hide()
# FIXME: There doesn't appear to be a good way to remove and restore
# toolbars wit having to rebuild them when they need to be displayed
# again.  So for now we'll leave the tab on the toolbox and simply
# hide the toolbar.  Which is a bit ugly, unfortunately.
# self.context.toolbox.remove_toolbar(2) # remove Where toolbar
      self.context.where_toolbar.hide()
#     self.context.toolbox.remove_toolbar(2) # remove When toolbar
      self.context.when_toolbar.hide()
    else:
      button4.show()
      label6.show()
      rb7.show()
      rb8.show()
      rb9.show()
      rb10.show()
      rb11.show()
      rb12.show()
      self.context.where_toolbar.show()
      self.context.when_toolbar.show()


  def never_visible(self, dec):
    if (latitude >= 0):
      return ((latitude + dec) < 0)
    else:
      return ((latitude + dec) > 0)


  def not_located(self, never):
    self.location.clear()
    self.plotchart()
    if (never):
      self.context.identifyobject.set_label(_('This object is always below the horizon at your location.'))
    else:
      if (self.magnifying):
        self.context.identifyobject.set_label(_('This object is currently not in the field of view.'))
      else:
        self.context.identifyobject.set_label(_('This object is currently below your horizon.'))

    
  def located(self, x, y):
    self.location.set(x, y)
    self.plotchart()


# If the left mouse button is pressed, identify the object at the
# cursor position.  If the right mouse button is pressed either
# Magnify (if the whole sky is shown) or revert to whole sky.

  def pressing(self, widget, event):
# Ignore any press that's outside the display circle.
    dx = (self.margin - 2 + self.xoffset + self.diameter / 2) - event.x
    dy = (self.margin - 2 + self.yoffset + self.diameter / 2) - event.y
    r = sqrt(dx * dx + dy * dy)
    if (r <= self.diameter / 2.0):
      if (event.button == 3):
        self.magnify(event.x, event.y)
      elif (event.button == 1):
        self.pmap.identify(event.x, event.y)
      else:
        pass

  def get_symbol_name(self, i):
    if (i == 2):
      return _('Moon')
    if (i == 7):
      return _('Sun')
    else:
      (name, wbar, e, a, I, O, L0, dL) = planets[i]
      return name

# ------------------------------------------------------------------------------
#
# Input Error Reporting methods

  def lon_error(self):
# FIXME: gtk.gdk.Display.beep()
    self.context.identifyobject.set_label(_('Input Error: The longitude was not understood.'))

  def lat_error(self):    
# FIXME: gtk.gdk.Display.beep()
    self.context.identifyobject.set_label(_('Input Error: The latitude was not understood.'))

  def time_error(self):
# FIXME: gtk.gdk.Display.beep()
    self.context.identifyobject.set_label(_('Input Error: The time was not understood.'))

  def zone_error(self):
# FIXME: gtk.gdk.Display.beep()
    self.context.identifyobject.set_label(_('Input Error: The time zone offset was not understood.'))


# -------------------------------------------------------------------------------
#
#   Methods for drawing the chart:

  def plotchart(self):
    if (self.canplot):
      self.plotfield()
      if (self.magnifying):
        self.plot_magnified()
      else:
        self.plot_whole_sky()
    return True


  def plotfield(self):
    global now
    global zoneoffset

# Erase prior plot

    if (not self.canplot):
      return
    self.cleararea()
    if (invertdisplay):
      if (nightvision):
        self.gc.set_foreground(self.colors[2])
      else:
        self.gc.set_foreground(self.colors[0])
    else:
      self.gc.set_foreground(self.colors[1])
    self.window.draw_arc(self.gc,
                              True,
                              self.xoffset + self.margin - 2,
                              self.yoffset + self.margin - 2,
                              self.diameter + 4,
                              self.diameter + 4,
                              0,
                              23040)

# Erase the pixels-to-object and object-to-pixels maps, since objects may now
# occupy different (x, y).

    self.pmap.clear()
    self.omap.clear()

# Plot sky circle

    if (not invertdisplay):
      if (nightvision):
        self.gc.set_foreground(self.colors[2])
      else:
        self.gc.set_foreground(self.colors[0])
    else:
      self.gc.set_foreground(self.colors[1])
    self.window.draw_arc(self.gc,
                              False,
                              self.xoffset + self.margin - 2,
                              self.yoffset + self.margin - 2,
                              self.diameter + 4,
                              self.diameter + 4,
                              0,
                              23040)

# label the cardinal points.

    if (nightvision):
      self.gc.set_foreground(self.colors[2])
    else:
      self.gc.set_foreground(self.colors[1])
    self.pangolayout.set_text(_('N'))
    self.window.draw_layout(self.gc,
                     self.xoffset + self.margin + self.diameter / 2 - 10,
                     self.margin - 30, self.pangolayout)
    self.pangolayout.set_text(_('S'))
    self.window.draw_layout(self.gc,
                     self.xoffset + self.margin + self.diameter / 2 - 10,
                     2 * self.margin + self.diameter - 30, self.pangolayout)
    if (not fliphorizontally):
      self.pangolayout.set_text(_('E'))
    else:
      self.pangolayout.set_text(_('W'))
    self.window.draw_layout(self.gc,
                     self.xoffset + self.margin - 30,
                     self.margin + self.diameter / 2 - 10, self.pangolayout)
    if (not fliphorizontally):
      self.pangolayout.set_text(_('W'))
    else:
      self.pangolayout.set_text(_('E'))
    self.window.draw_layout(self.gc,
                     self.xoffset + self.margin + self.diameter + 10,
                     self.margin + self.diameter / 2 - 10, self.pangolayout)
    if (not invertdisplay):
      if (nightvision):
        self.gc.set_foreground(self.colors[2])
      else:
        self.gc.set_foreground(self.colors[0])
    else:
      self.gc.set_foreground(self.colors[1])

# Set the time of plotting (now).

    if (not specifytime):
      now = datetime.utcnow()
      (tstr, ostr) = set_time_and_UTC_offset()
      entry3.set_text(tstr)
      entry4.set_text(ostr)
    else:
      now = get_time_and_UTC_offset(entry3.get_text(), entry4.get_text())
      (hh, mm) = parse_zone_offset(entry4.get_text())
      zoneoffset = 60 * hh
      if (hh < 0):
        zoneoffset = zoneoffset - mm
      else:
        zoneoffset = zoneoffset + mm
    return True


  def plot_whole_sky(self):
      
# Plot the entire visible sky, starting with the stars.

    self.plot_all_stars()
    self.plot_all_DSOs()
    self.plot_all_constellations()
    self.plot_all_planets()
    return True


  def plot_magnified(self):
# Plot a circular section of the sky 3.5 degrees in radius centered about the
# coordinates saved in self.mag_center.
    self.plot_mag_stars()
    self.plot_mag_DSOs()
# we don't plot constellation figures in the magnified view.
    self.plot_mag_planets()
    return True


  def plot_all_stars(self):
    for name, (ra, dec, mag, cid) in star_chart.iteritems():

# convert the ra and dec from the J2000 epoch to the plot time

      polar = epochpolartonow((ra, dec), now)

# convert the equatorial coordinates to altitude and azimuth

      azalt = polartoazalt(polar, now)

# plot any object that is more than 1 degree above the horizon

      if (azalt[1] >= 0.0175):
        starsize = 2 + int(7.0 - mag)
        (px, py) = self.azalttoxy(azalt)
        px = px + self.margin - 2 + self.xoffset - starsize / 2
        py = py + self.margin - 2 + self.yoffset - starsize / 2

# Add the star to the omap even if it was too dim to plot.

        self.omap.add('star', name, px, py)

# if the star is bright enough, add it to pmap and plot it.

        if (mag <= limitingmagnitude):
          self.pmap.add(px, py, 'star', name)
          self.plot_star(px, py, starsize)


  def plot_all_DSOs(self):

# Plot the deep sky objects.

    for i in range(len(dso_chart)):
      (nM, strCon, ra, dec, mag, majA, minA, posA, strT, strN) = dso_chart[i]

# convert the ra and dec from the J2000 epoch to the plot time

      polar = epochpolartonow((ra, dec), now)

# convert the equatorial coordinates to altitude and azimuth

      azalt = polartoazalt(polar, now)
      if (azalt[1] >= 0.0175):
        (px, py) = self.azalttoxy(azalt)
        px = px + self.margin - 2 + self.xoffset
        py = py + self.margin - 2 + self.yoffset
        self.plot_DSO(strT, majA, minA, mag, px, py)

# Add the DSO to the maps.

        if (strN == ''):
          if (nM[0:1] == 'M'):
            self.pmap.add(px, py, 'dso-messier', nM)
            self.omap.add('dso-messier', nM, px, py)
          elif (nM[0:1] == 'N'):
            self.pmap.add(px, py, 'dso-ngc', nM)
            self.omap.add('dso-ngc', nM, px, py)
          else:
            self.pmap.add(px, py, 'dso', nM)
            self.omap.add('dso', nM, px, py)
        else:
          if (nM[0:1] == 'M'):
            self.pmap.add(px, py, 'dso-messier', strN + ' (' + nM +')')
            self.omap.add('dso-messier',  strN + ' (' + nM + ')', px, py)
          elif (nM[0:1] == 'N'):
            self.pmap.add(px, py, 'dso-ngc', strN + ' (' + nM +')')
            self.omap.add('dso-ngc',  strN + ' (' + nM + ')', px, py)
          else:
            self.pmap.add(px, py, 'dso', strN + ' (' + nM +')')
            self.omap.add('dso',  strN + ' (' + nM + ')', px, py)


  def plot_all_constellations(self):

# Plot the constellation figures.  This is essentially the same process as for
# plotting a star but we have to figure out the alt/az coordinates for both ends
# of the line segment.

    if (drawconstellations):
      if (not invertdisplay):
        if (nightvision):
          self.gc.set_foreground(self.colors[2])
        else:
          self.gc.set_foreground(self.colors[0])
      else:
        self.gc.set_foreground(self.colors[1])
      for code, (name, lines) in figures.iteritems():
        for i in range(len(lines)):
          (ra1, dec1, ra2, dec2) = lines[i]
          polar1 = epochpolartonow((ra1, dec1), now)
          azalt1 = polartoazalt(polar1, now)
          if (azalt1[1] >=  0.0175):
            (px1, py1) = self.azalttoxy(azalt1)
            px1 = px1 + self.margin - 2 + self.xoffset
            py1 = py1 + self.margin - 2 + self.yoffset
            polar2 = epochpolartonow((ra2, dec2), now)
            azalt2 = polartoazalt(polar2, now)
            if (azalt2[1] >=  0.0175):
              (px2, py2) = self.azalttoxy(azalt2)
              px2 = px2 + self.margin - 2 + self.xoffset
              py2 = py2 + self.margin - 2 + self.yoffset
              self.window.draw_line(self.gc, px1, py1, px2, py2)


  def plot_all_planets(self):

# Plot the planets, the moon and the sun.

# We need to adjust the Keplerian data by the number of centuries since the
# epoch.  Mostly, this causes the planets to revolve around the sun at the rate
# of their orbital period (N in the following computations), but the orbits
# themselves will also change over time.
# FIXME: for now, we only worry about orbital motion, not additionally orbital
#        mutation.

    T = (jtime(now) - 2451545.0) / 36525.0

# Calculate the earth's heliocentric radius and longitude coordinates so we
# can figure the geocentric ecliptic coordinates of the other planets.

    (name, wbar, e, a, I, O, L0, dL) = planets[2]
    N = dL * T
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0
# FIXME: Should adjust the orbit's inclination and orientation
#    I = I0 + dI * T
#    wbar = wbar0 + dwbar * T
#    O = O0 + dO * T
    M = N + L0 - wbar
    while (M >= 360.0):
      M = M - 360.0
    while (M < 0.0):
      M = M + 360.0
    Le = N + 2 * rtod(e) * sin(dtor(M)) + L0
    while (Le >= 360.0):
      Le = Le - 360.0
    while (Le < 0.0):
      Le = Le + 360.0
    v = Le - wbar
    Re = a * (1 - e * e) / (1 + e * cos(dtor(v)))

# now plot the planets.

    for i in range(len(planets)):
      if (i == 2):
        continue
      (name, wbar, e, a, I, O, L0, dL) = planets[i]
      N = dL * T
      while (N >= 360.0):
        N = N - 360.0
      while (N < 0.0):
        N = N + 360.0
# FIXME:  For Mercury, especially, dI, dwbar and dO are large enough to matter.
#      I = I0 + dI * T
#      wbar = wbar0 + dwbar * T
#      O = O0 + dO * T

# compute the mean anomaly

      M = N + L0 - wbar
      while (M >= 360.0):
        M = M - 360.0
      while (M < 0.0):
        M = M + 360.0

# compute the heliocentric longitude
# FIXME: Rather than solve Keppler's equation (M = E - 180/pi * e * sin(E)) for
#        E, we will use the approximation E ~ M and calculate the heliocentric
#        longitude using the mean anomaly instead of the eccentric anomaly.
#        This approximation is pretty close except for Mercury and Pluto which
#        have values of e high enough to make a difference.

      l = N + 2 * rtod(e) * sin(dtor(M)) + L0
      while (l >= 360.0):
        l = l - 360.0
      while (l < 0.0):
        l = l + 360.0

# now calculate the actual anomaly

      v = l - wbar

# find the planet's radial distance

      r = a * (1 - e * e) / (1 + e * cos(dtor(v)))

# calculate the heliocentric latitude

      psi = rtod(asin(sin(dtor(l - O)) * sin(dtor(I))))

# project to the ecliptic

      y = sin(dtor(l - O)) * cos(dtor(I))
      x = cos(dtor(l - O))
      k = rtod(atan2(y, x))
      lprime = k + O
      while (lprime >= 360.0):
        lprime = lprime - 360.0
      while (lprime < 0.0):
        lprime = lprime + 360.0
      rprime = r * cos(dtor(psi))

# Using the coordinates we already calculated for the current position of the
# earth, convert the above coordinates to geocentric latitude and longitude

      if (i < 2):

# for inner planets, use this equation to get longitude:

        y = rprime * sin(dtor(Le - lprime))
        x = Re - rprime * cos(dtor(Le - lprime))
        k = rtod(atan2(y, x))
        lam = 180.0 + Le + k 
      else:

# for outer planets, use this equation to get longitude:

        y = Re * sin(dtor(lprime - Le))
        x = rprime - Re * cos(dtor(lprime - Le))
        k = rtod(atan2(y, x))
        lam = k + lprime
      while (lam >= 360.0):
        lam = lam - 360.0
      while (lam < 0.0):
        lam = lam + 360.0

# all planets use the same equation for ecliptic latitude (and this atan has no quadrant ambiguity):

      y = rprime * tan(dtor(psi)) * sin(dtor(lam - lprime))
      x = Re * sin(dtor(lprime - Le))
      beta = rtod(atan(y/x))

# convert lam and beta to RA and DEC

      y = sin(dtor(lam)) * cos(eps) - tan(dtor(beta)) * sin(eps)
      x = cos(dtor(lam))
      k = rtod(atan2(y, x))
      ra = k / 15.0
      while (ra < 0.0):
        ra = ra + 24.0
      dec = rtod(asin(sin(dtor(beta)) * cos(eps) + cos(dtor(beta)) * \
                      sin(eps) * sin(dtor(lam))))
# remember this planet's declination.  Locate needs to know.
      dec_from_planet[name] = dec
      
# convert to azalt

      azalt = polartoazalt((ra, dec), now)

# convert to x,y and plot if the planet is above the horizon

      if (azalt[1] >= 0.0175):
        (px, py) = self.azalttoxy(azalt)
        px = px + self.margin - 2 + self.xoffset
        py = py + self.margin - 2 + self.yoffset
        self.plot_planetary_symbol(i, px, py)

# Add the planet to the maps

        self.pmap.add(px, py, 'planet', name)
        self.omap.add('planet', name, px, py)

# Plot the sun.  This is virtually the same as for a planet, but we can simplify
# some computations because the sun is by definition on the ecliptic.

    (name, wbar, e, a, I, O, L0, dL) = sun
    N = dL * T
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0    
    M = N + L0 - wbar
    while (M >= 360.0):
      M = M - 360.0
    while (M < 0.0):
      M = M + 360.0
    v = M + 2 * rtod(e) * sin(dtor(M))
    while (v >= 360.0):
      v = v - 360.0
    while (lam < 0.0):
      v = v + 360.0
    lam = v + wbar
    while (lam >= 360.0):
      lam = lam - 360.0
    while (lam < 0.0):
      lam = lam + 360.0
    y = sin(dtor(lam)) * cos(eps)
    x = cos(dtor(lam))
    k = rtod(atan2(y, x))
    ra = k / 15.0
    while (ra < 0.0):
      ra = ra + 24.0

# because beta is (by definition) 0, calculating dec is far simpler:

    dec = rtod(asin(sin(sin(eps) * sin(dtor(lam)))))
# remember this planet's declination.  Locate needs to know.
    dec_from_planet[name] = dec
    azalt = polartoazalt((ra, dec), now)
    if (azalt[1] >= 0.0175):
      (px, py) = self.azalttoxy(azalt)
      px = px + self.margin - 2 + self.xoffset
      py = py + self.margin - 2 + self.yoffset
      self.plot_planetary_symbol(7, px, py)

# Add the sun to maps

      self.pmap.add(px, py, 'planet', name)
      self.omap.add('planet', name, px, py)

# Plot the moon.  Since the moon's orbit varies radically over time, there are
# a lot of 'fudge factor' correction terms in these computations.

    (name, L0, P0, N0, I, e, a, phi0, tau) = moon
    D = jtime(now) - float(tau)
    N = 360.0 / 365.242191 * D
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0
    M = N + L0 - P0
    while (M >= 360.0):
      M = M - 360.0
    while (M < 0.0):
      M = M + 360.0
    Ec = 2 * rtod(e) * sin(dtor(M))
    lam = Ec + L0
    while (lam >= 360.0):
      lam = lam - 360.0
    while (lam < 0.0):
      lam = lam + 360.0
    l = 13.176396 * D + L0               
    while (l >= 360.0):
      l = l - 360.0
    while (l < 0.0):
      l = l + 360.0
    Mm = l - 0.1114041 * D - P0
    while (Mm >= 360.0):
      Mm = Mm - 360.0
    while (Mm < 0.0):
      Mm = Mm + 360.0
    N = N0 - 0.0529539 * D
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0
    Ev = 1.2739 * sin(dtor(2 * (1 - lam) - Mm))
    Ae = 0.1858 * sin(dtor(M))
    A3 = 0.37 * sin(dtor(M))
    Mprime = Mm + Ev - Ae - A3
    Ec = 6.2886 * sin(dtor(Mprime))
    A4 = 0.214 * sin (dtor(2 * Mprime))
    lprime = l + Ev  + Ec - Ae + A4
    V = 0.6583 * sin(dtor(2 * (lprime - lam)))
    lon = lprime + V
    Nprime = N - 0.16 * sin(dtor(M))
    y = sin(dtor(lon - Nprime)) * cos(dtor(I))
    x = cos(dtor(lon - Nprime))
    k = rtod(atan2(y, x))
    lamM = Nprime + k
    while (lamM >= 360.0):
      lamM = lamM - 360.0
    while (lamM < 0.0):
      lamM = lamM + 360.0
    betM = rtod(asin(sin(dtor(lon - Nprime)) * sin(dtor(I))))
    y = sin(dtor(lamM)) * cos(eps) - tan(dtor(betM)) * sin(eps)
    x = cos(dtor(lamM))
    k = rtod(atan2(y, x))
    ra = k / 15.0
    while (ra < 0.0):
      ra = ra + 24.0
    while (ra > 24.0):
      ra = ra - 24.0
    dec = rtod(asin(sin(dtor(betM)) * cos(eps) + cos(dtor(betM)) * \
                    sin(eps) * sin(dtor(lamM))))
# FIXME: We aren't accounting for parallax: RA and DEC are actual, not
#	 apparent.  This could introduce as much as a degree (approx.)
#	 of error.  Not a big deal if one is not using a telescope to
#	 observe the moon, except when the moon's position relative to
#	 a planet or star (as in conjunctions and occultations)
#	 matters.

# remember this planet's declination.  Locate needs to know.
    dec_from_planet[name] = dec

# convert to azalt

    azalt = polartoazalt((ra, dec), now)
    if (azalt[1] >= 0.0175):
      (px, py) = self.azalttoxy(azalt)
      px = px + self.margin - 2 + self.xoffset
      py = py + self.margin - 2 + self.yoffset
      self.plot_planetary_symbol(2, px, py)

# Add the moon to maps

      self.pmap.add(px, py, 'planet', name)
      self.omap.add('planet', name, px, py)
    self.gc.set_foreground(self.colors[1])
    self.location.plot_cross()
    return True


  def vector_from_fov(self, ra, dec, now):

# The center of the FOV is self.mag_center[]; [0] is its x; [1] is its y.  These
# coordinates are screen-absolute.  Make them center-of-plot-relative.

    x0 = self.mag_center[0] - (self.margin - 2 + self.xoffset + self.diameter / 2)
#    if (fliphorizontally):
#      x0 = (self.margin - 2 + self.xoffset + self.diameter / 2) - self.mag_center[0]
#    else:
#      x0 = self.mag_center[0] - (self.margin - 2 + self.xoffset + self.diameter / 2)
    y0 = (self.margin - 2 + self.yoffset + self.diameter / 2) - self.mag_center[1]
    y0 = self.mag_center[1] - (self.margin - 2 + self.yoffset + self.diameter / 2)

# Convert ra, dec and now to az and alt and then to x and y.

    (az, alt) = polartoazalt((ra, dec), now)
    if (alt < 0.0):
      return (-9999, -9999) # object cannot possibly be seen.
    (x, y) = self.azalttoxy((az, alt))

# Determine the difference between mag_center[] and (x,y).

    x = x - self.diameter / 2
    y = y - self.diameter / 2
    dx = x - x0
    dy = y - y0

# Apply magnification factor (about 25X; the field of 7x binoculars is around 7 degrees).

    dx = dx * 25
    dy = dy * 25

# Is this object within the plot circle?

    r = sqrt(dx * dx + dy * dy)
    if (r > self.diameter / 2.0):
      return (-9999, -9999) # object cannot possibly be seen.

# Return dx,dy -- this is the coordinate of the object relative to the
# plot center.

    return (dx, dy)


  def plot_mag_stars(self):

# Plot the stars within the field of view

    for name, (ra, dec, mag, cid) in star_chart.iteritems():

# Convert the ra and dec from the J2000 epoch to the plot time

      (ra1, dec1) = epochpolartonow((ra, dec), now)
      starsize = 2 + int(9.0 - mag)
      (px, py) = self.vector_from_fov(ra1, dec1, now)
      if (px < -9000) or (py < -9000):
        continue # object is outside of the FOV or within 0.5 percent of it.
      px = int(px) + self.margin - 2 + self.xoffset + self.diameter / 2 - starsize / 2
      py = int(py) + self.margin - 2 + self.yoffset + self.diameter / 2 - starsize / 2

# Add the star to the omap even if it was too dim to plot.

      self.omap.add('star', name, px, py)

# All stars are bright enough in magnified mode, so add this star to
# pmap and plot it.

      self.pmap.add(px, py, 'star', name)
      self.plot_star(px, py, starsize)


  def plot_mag_DSOs(self):
 
# Plot only those DSOs that are within 3.5 degrees of the
# magnification coordinates.

    for i in range(len(dso_chart)):
      (nM, strCon, ra, dec, mag, majA, minA, posA, strT, strN) = dso_chart[i]

# Convert the ra and dec from the J2000 epoch to the plot time

      (ra1, dec1) = epochpolartonow((ra, dec), now)
      (px, py) = self.vector_from_fov(ra1, dec1, now)
      if (px < -9000) or (py < -9000):
        continue # object is outside of the FOV or within 0.5 percent of it.
      px = int(px) + self.margin - 2 + self.xoffset + self.diameter / 2
      py = int(py) + self.margin - 2 + self.yoffset + self.diameter / 2
      self.plot_DSO(strT, majA, minA, mag, px, py)

# Add the DSO to the maps.

      if (strN == ''):
        if (nM[0:1] == 'M'):
          self.pmap.add(px, py, 'dso-messier', nM)
          self.omap.add('dso-messier', nM, px, py)
        elif (nM[0:1] == 'N'):
          self.pmap.add(px, py, 'dso-ngc', nM)
          self.omap.add('dso-ngc', nM, px, py)
        else:
          self.pmap.add(px, py, 'dso', nM)
          self.omap.add('dso', nM, px, py)
      else:
        if (nM[0:1] == 'M'):
          self.pmap.add(px, py, 'dso-messier', strN + ' (' + nM +')')
          self.omap.add('dso-messier',  strN + ' (' + nM + ')', px, py)
        elif (nM[0:1] == 'N'):
          self.pmap.add(px, py, 'dso-ngc', strN + ' (' + nM +')')
          self.omap.add('dso-ngc',  strN + ' (' + nM + ')', px, py)
        else:
          self.pmap.add(px, py, 'dso', strN + ' (' + nM +')')
          self.omap.add('dso',  strN + ' (' + nM + ')', px, py)


  def plot_mag_planets(self):

# Plot the planets which are within the field of view

    T = (jtime(now) - 2451545.0) / 36525.0

# Calculate the earth's heliocentric radius and longitude coordinates so we
# can figure the geocentric ecliptic coordinates of the other planets.

    (name, wbar, e, a, I, O, L0, dL) = planets[2]
    N = dL * T
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0
# FIXME: Should adjust the orbit's inclination and orientation
#    I = I0 + dI * T
#    wbar = wbar0 + dwbar * T
#    O = O0 + dO * T
    M = N + L0 - wbar
    while (M >= 360.0):
      M = M - 360.0
    while (M < 0.0):
      M = M + 360.0
    Le = N + 2 * rtod(e) * sin(dtor(M)) + L0
    while (Le >= 360.0):
      Le = Le - 360.0
    while (Le < 0.0):
      Le = Le + 360.0
    v = Le - wbar
    Re = a * (1 - e * e) / (1 + e * cos(dtor(v)))

# now plot the planets.

    for i in range(len(planets)):
      if (i == 2):
        continue
      (name, wbar, e, a, I, O, L0, dL) = planets[i]
      N = dL * T
      while (N >= 360.0):
        N = N - 360.0
      while (N < 0.0):
        N = N + 360.0
# FIXME:  For Mercury, especially, dI, dwbar and dO are large enough to matter.
#      I = I0 + dI * T
#      wbar = wbar0 + dwbar * T
#      O = O0 + dO * T

# compute the mean anomaly

      M = N + L0 - wbar
      while (M >= 360.0):
        M = M - 360.0
      while (M < 0.0):
        M = M + 360.0

# compute the heliocentric longitude
# FIXME: Rather than solve Keppler's equation (M = E - 180/pi * e * sin(E)) for
#        E, we will use the approximation E ~ M and calculate the heliocentric
#        longitude using the mean anomaly instead of the eccentric anomaly.
#        This approximation is pretty close except for Mercury and Pluto which
#        have values of e high enough to make a difference.

      l = N + 2 * rtod(e) * sin(dtor(M)) + L0
      while (l >= 360.0):
        l = l - 360.0
      while (l < 0.0):
        l = l + 360.0

# now calculate the actual anomaly

      v = l - wbar

# find the planet's radial distance

      r = a * (1 - e * e) / (1 + e * cos(dtor(v)))

# calculate the heliocentric latitude

      psi = rtod(asin(sin(dtor(l - O)) * sin(dtor(I))))

# project to the ecliptic

      y = sin(dtor(l - O)) * cos(dtor(I))
      x = cos(dtor(l - O))
      k = rtod(atan2(y, x))
      lprime = k + O
      while (lprime >= 360.0):
        lprime = lprime - 360.0
      while (lprime < 0.0):
        lprime = lprime + 360.0
      rprime = r * cos(dtor(psi))

# Using the coordinates we already calculated for the current position of the
# earth, convert the above coordinates to geocentric latitude and longitude

      if (i < 2):

# for inner planets, use this equation to get longitude:

        y = rprime * sin(dtor(Le - lprime))
        x = Re - rprime * cos(dtor(Le - lprime))
        k = rtod(atan2(y, x))
        lam = 180.0 + Le + k 
      else:

# for outer planets, use this equation to get longitude:

        y = Re * sin(dtor(lprime - Le))
        x = rprime - Re * cos(dtor(lprime - Le))
        k = rtod(atan2(y, x))
        lam = k + lprime
      while (lam >= 360.0):
        lam = lam - 360.0
      while (lam < 0.0):
        lam = lam + 360.0

# All planets use the same equation for ecliptic latitude (and this
# atan has no quadrant ambiguity):

      y = rprime * tan(dtor(psi)) * sin(dtor(lam - lprime))
      x = Re * sin(dtor(lprime - Le))
      beta = rtod(atan(y/x))

# Convert lam and beta to RA and DEC

      y = sin(dtor(lam)) * cos(eps) - tan(dtor(beta)) * sin(eps)
      x = cos(dtor(lam))
      k = rtod(atan2(y, x))
      ra1 = k / 15.0
      while (ra1 < 0.0):
        ra1 = ra1 + 24.0
      dec1 = rtod(asin(sin(dtor(beta)) * cos(eps) + cos(dtor(beta)) * \
                      sin(eps) * sin(dtor(lam))))

# Remember this planet's declination.  Locate needs to know.

      dec_from_planet[name] = dec1
      (px, py) = self.vector_from_fov(ra1, dec1, now)
      if (px < -9000) or (py < -9000):
        continue # object is outside of the FOV or within 0.5 percent of it.
      px = int(px) + self.margin - 2 + self.xoffset + self.diameter / 2
      py = int(py) + self.margin - 2 + self.yoffset + self.diameter / 2
      self.plot_planetary_symbol(i, px, py)

# Add the planet to the maps

      self.pmap.add(px, py, 'planet', name)
      self.omap.add('planet', name, px, py)

# Plot the sun.  This is virtually the same as for a planet, but we can simplify
# some computations because the sun is by definition on the ecliptic.

    (name, wbar, e, a, I, O, L0, dL) = sun
    N = dL * T
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0    
    M = N + L0 - wbar
    while (M >= 360.0):
      M = M - 360.0
    while (M < 0.0):
      M = M + 360.0
    v = M + 2 * rtod(e) * sin(dtor(M))
    while (v >= 360.0):
      v = v - 360.0
    while (lam < 0.0):
      v = v + 360.0
    lam = v + wbar
    while (lam >= 360.0):
      lam = lam - 360.0
    while (lam < 0.0):
      lam = lam + 360.0
    y = sin(dtor(lam)) * cos(eps)
    x = cos(dtor(lam))
    k = rtod(atan2(y, x))
    ra1 = k / 15.0
    while (ra1 < 0.0):
      ra1 = ra1 + 24.0

# Because beta is (by definition) 0, calculating dec is far simpler:

    dec1 = rtod(asin(sin(sin(eps) * sin(dtor(lam)))))
    (px, py) = self.vector_from_fov(ra1, dec1, now)
    if (px >= -9000) and (py >= -9000):
# Object is not outside of the FOV.
      px = int(px) + self.margin - 2 + self.xoffset + self.diameter / 2
      py = int(py) + self.margin - 2 + self.yoffset + self.diameter / 2
      self.plot_planetary_symbol(7, px, py)

# Add the sun to maps

      self.pmap.add(px, py, 'planet', name)
      self.omap.add('planet', name, px, py)

# Plot the moon.  Since the moon's orbit varies radically over time, there are
# a lot of "fudge factor" correction terms in these computations.

    (name, L0, P0, N0, I, e, a, phi0, tau) = moon
    D = jtime(now) - float(tau)
    N = 360.0 / 365.242191 * D
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0
    M = N + L0 - P0
    while (M >= 360.0):
      M = M - 360.0
    while (M < 0.0):
      M = M + 360.0
    Ec = 2 * rtod(e) * sin(dtor(M))
    lam = Ec + L0
    while (lam >= 360.0):
      lam = lam - 360.0
    while (lam < 0.0):
      lam = lam + 360.0
    l = 13.176396 * D + L0               
    while (l >= 360.0):
      l = l - 360.0
    while (l < 0.0):
      l = l + 360.0
    Mm = l - 0.1114041 * D - P0
    while (Mm >= 360.0):
      Mm = Mm - 360.0
    while (Mm < 0.0):
      Mm = Mm + 360.0
    N = N0 - 0.0529539 * D
    while (N >= 360.0):
      N = N - 360.0
    while (N < 0.0):
      N = N + 360.0
    Ev = 1.2739 * sin(dtor(2 * (1 - lam) - Mm))
    Ae = 0.1858 * sin(dtor(M))
    A3 = 0.37 * sin(dtor(M))
    Mprime = Mm + Ev - Ae - A3
    Ec = 6.2886 * sin(dtor(Mprime))
    A4 = 0.214 * sin (dtor(2 * Mprime))
    lprime = l + Ev  + Ec - Ae + A4
    V = 0.6583 * sin(dtor(2 * (lprime - lam)))
    lon = lprime + V
    Nprime = N - 0.16 * sin(dtor(M))
    y = sin(dtor(lon - Nprime)) * cos(dtor(I))
    x = cos(dtor(lon - Nprime))
    k = rtod(atan2(y, x))
    lamM = Nprime + k
    while (lamM >= 360.0):
      lamM = lamM - 360.0
    while (lamM < 0.0):
      lamM = lamM + 360.0
    betM = rtod(asin(sin(dtor(lon - Nprime)) * sin(dtor(I))))
    y = sin(dtor(lamM)) * cos(eps) - tan(dtor(betM)) * sin(eps)
    x = cos(dtor(lamM))
    k = rtod(atan2(y, x))
    ra1 = k / 15.0
    while (ra1 < 0.0):
      ra1 = ra1 + 24.0
    while (ra1 > 24.0):
      ra1 = ra1 - 24.0
    dec1 = rtod(asin(sin(dtor(betM)) * cos(eps) + cos(dtor(betM)) * \
                    sin(eps) * sin(dtor(lamM))))
# FIXME: We aren't accounting for parallax: RA and DEC are actual, not apparent.
#	 This could introduce as much as a degree (approx.) of error.  Not a big deal if
#        one is not using a telescope to observe the moon, except when the
#        moon's position relative to a planet or star (as in conjunctions and
#        occultations) matters.

# remember this planet's declination.  Locate needs to know.
    dec_from_planet[name] = dec1
    (px, py) = self.vector_from_fov(ra1, dec1, now)
    if (px >= -9000) and (py >= -9000):
# object is not outside of the FOV.
      px = int(px) + self.margin - 2 + self.xoffset + self.diameter / 2
      py = int(py) + self.margin - 2 + self.yoffset + self.diameter / 2
      self.plot_planetary_symbol(2, px, py)

# Add the moon to maps

      self.pmap.add(px, py, 'planet', name)
      self.omap.add('planet', name, px, py)
    self.gc.set_foreground(self.colors[1])
    self.location.plot_cross()
    return True


  def plot_star(self, px, py, starsize):
    self.window.draw_arc(self.gc, True,
             px,
             py,
             starsize,
             starsize,
             0,
             360*64)


  def plot_planetary_symbol(self, i, px, py):

# i is planet number (0 = mercury, 1= venus, 2 = moon, 3 = mars, 4 = jupiter,
# 5 = saturn, 6 = uranus, 7 = sun); (px, py) is the center point of the symbol.

    if (i == 0):

# mercury

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, False, px-5, py-7, 10, 10, 0, 360*64)
      self.window.draw_line(self.gc, px+4, py-9, px+4, py-7)
      self.window.draw_line(self.gc, px-4, py-9, px-4, py-7)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_line(self.gc, px, py+3, px, py+7)
      self.window.draw_line(self.gc, px-2, py+5, px+2, py+5)
    elif (i == 1):

# venus

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, False, px-5, py-7, 10, 10, 0, 360*64)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_line(self.gc, px, py+3, px, py+7)
      self.window.draw_line(self.gc, px-2, py+5, px+2, py+5)
    elif (i == 2):

# moon

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_polygon(self.gc, True, ((px+1, py-11), (px+4, py-11),
                                               (px+5, py-10), (px+6, py-9),
                                               (px+7, py-8),  (px+8, py-7),
                                               (px+10, py-2), (px+12,py),
                                               (px+10, py+2), (px+8, py+7),
                                               (px+7, py+8),  (px+6, py+9),
                                               (px+5, py+10), (px+4,py+11),
                                               (px+1, py+11), (px+4, py+4),
                                               (px+6, py+2),  (px+6, py-2),
                                               (px+4, py-4)))
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
    elif (i == 3):

# mars

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, False, px-6, py-4, 10, 10, 0, 360*64)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_line(self.gc, px+2, py-2, px+6, py-6)
      self.window.draw_line(self.gc, px+3, py-6, px+6, py-6)
      self.window.draw_line(self.gc, px+6, py-6, px+6, py-3)
    elif (i == 4):

# jupiter

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_line(self.gc, px-6, py-6, px-4, py-8)
      self.window.draw_line(self.gc, px-4, py-8, px-2, py-8)
      self.window.draw_line(self.gc, px-2, py-8, px+1, py-6)
      self.window.draw_line(self.gc, px+1, py-6, px-5, py+2)
      self.window.draw_line(self.gc, px-5, py+2, px+7, py+2)
      self.window.draw_line(self.gc, px+4, py-8, px+4, py+7)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
    elif (i == 5):

# saturn

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_line(self.gc, px-6, py-6, px-6, py+5)
      self.window.draw_line(self.gc, px-6, py, px-5, py-1)
      self.window.draw_line(self.gc, px-5, py-1, px-4, py-2)
      self.window.draw_line(self.gc, px-4, py-2, px-1, py-3)
      self.window.draw_line(self.gc, px-1, py-3, px, py-4)
      self.window.draw_line(self.gc, px, py-4, px+1, py+1)
      self.window.draw_line(self.gc, px+1, py+1, px-1, py+4)
      self.window.draw_line(self.gc, px-1, py+4, px, py+5)
      self.window.draw_line(self.gc, px, py+5, px+6, py+4)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
    elif (i == 6):

# uranus

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, False, px-5, py-3, 10, 10, 0, 360*64)
      self.window.draw_arc(self.gc, True, px-2, py, 4, 4, 0, 360*64)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_line(self.gc, px, py-3, px, py-9)
      self.window.draw_line(self.gc, px-2, py-5, px, py-9)
      self.window.draw_line(self.gc, px+2, py-5, px, py-9)
    else:

# sun

      self.gc.set_line_attributes(2, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)
      self.window.draw_arc(self.gc, False, px-12, py-12, 24, 24, 0, 360*64)
      self.window.draw_arc(self.gc, True, px-2, py-2, 4, 4, 0, 360*64)
      self.gc.set_line_attributes(1, gtk.gdk.LINE_SOLID, gtk.gdk.CAP_BUTT,
                                  gtk.gdk.JOIN_MITER)


  def plot_DSO(self, type, maja, mina, mag, px, py):
    if (not invertdisplay):
      if (nightvision):
        fg_color = self.colors[2]
      else:
        fg_color = self.colors[0]
    else:
      fg_color = self.colors[1]
    if (mag > 9.0) or (maja < 10.0):
# this object is too dim or too small to be interesting
      pass
    else:
# FIXME:  The Messier Catalog should include the object's position angle and
#         this code should plot it in the proper orientation.  What we have here
#         will always plot the object with the major axis east-west and the
#         minor axis north-south.

# We plot these objects actual size.  Of course, they are mostly too small to plot.

      if (self.magnifying):
        dx = int(maja / 60.0 * self.diameter / 7.0)
        dy = int(mina / 60.0 * self.diameter / 7.0)
        if (dx < 2) or (dy < 2):
          return # too small.
      else:
        dx = int(maja / 60.0 * self.diameter / 180.0)
        dy = int(mina / 60.0 * self.diameter / 180.0)
        if (dx < 2) or (dy < 2):
          return # too small.
      if (type == 'Gal'):
# plot as gray ellipse with solid outline.
        self.gc.set_foreground(self.colors[3])
        self.window.draw_arc(self.gc,
                              True,
                              px - int(dx / 2),
                              py - int(dy / 2),
                              dx - 1,
                              dy - 1,
                              0,
                              23040)
        self.gc.set_foreground(fg_color)
        self.window.draw_arc(self.gc,
                              False,
                              px - int(dx / 2),
                              py - int(dy / 2),
                              dx,
                              dy,
                              0,
                              23040)
      elif (type == 'PlN'):
# plot as gray circle with central dot
        self.gc.set_foreground(self.colors[3])
        self.window.draw_arc(self.gc,
                              True,
                              px - int(dx / 2),
                              py - int(dx / 2),
                              dx,
                              dx,
                              0,
                              23040)
        self.gc.set_foreground(fg_color)
        self.window.draw_arc(self.gc,
                              True,
                              px - 2,
                              py - 2,
                              4,
                              4,
                              0,
                              23040)
      elif (type == 'SNR') or (type == 'OCl'):
# plot as gray circle with no outline.
        self.gc.set_foreground(self.colors[3])
        self.window.draw_arc(self.gc,
                              True,
                              px - int(dx / 2),
                              py - int(dx / 2),
                              dx,
                              dx,
                              0,
                              23040)
        self.gc.set_foreground(fg_color)
      elif (type == 'C/N') or (type == 'DfN'):
# plot as gray rectangle with no outline.
        self.gc.set_foreground(self.colors[3])
        self.window.draw_rectangle(self.gc,
                                        True,
                                        px - int(dx / 2),
                                        py - int(dy / 2),
                                        dx,
                                        dy)
        self.gc.set_foreground(fg_color)
      elif (type == 'GCl'):
# plot as gray circle with outline and central dot.
        self.gc.set_foreground(self.colors[3])
        self.window.draw_arc(self.gc,
                              True,
                              px - int(dx / 2),
                              py - int(dx / 2),
                              dx,
                              dx,
                              0,
                              23040)
        self.gc.set_foreground(fg_color)
        self.window.draw_arc(self.gc,
                              False,
                              px - int(dx / 2),
                              py - int(dx / 2),
                              dx - 1,
                              dx - 1,
                              0,
                              23040)
        self.window.draw_arc(self.gc,
                              True,
                              px - 2,
                              py - 2,
                              4,
                              4,
                              0,
                              23040)
      else:
#	Dbl = double star
#	??? = unknown or unclassified object
# these are not plotted.
        pass

      
  def cleararea(self):
    
# Clear the drawing surface

    if (nightvision):
      self.gc.set_foreground(self.colors[1])
    else:
      self.gc.set_foreground(self.colors[3])
    self.window.draw_rectangle(self.gc,
                                    True,
                                    1,
                                    1,
                                    self.screensize[0],
                                    self.screensize[1])
    label1.queue_draw()
    label2.queue_draw()
    label4.queue_draw()
    label5.queue_draw()
    label6.queue_draw()
    button1.queue_draw()
    button2.queue_draw()
    button3.queue_draw()
    button4.queue_draw()
    button5.queue_draw()
    button51.queue_draw()
    button6.queue_draw()
    rb1.queue_draw()
    rb2.queue_draw()
    rb3.queue_draw()
    rb4.queue_draw()
    rb5.queue_draw()
    rb6.queue_draw()
    rb7.queue_draw()
    rb8.queue_draw()
    rb9.queue_draw()
    rb10.queue_draw()
    rb11.queue_draw()
    rb12.queue_draw()


# ============================== StarChart Object =============================

class StarChart(activity.Activity):
  def __init__(self, handle):
    global now
    global zoneoffset
    global abbrev_from_name
    global longitude
    global latitude
    activity.Activity.__init__(self, handle)
    os.chdir(get_bundle_path())
    self.set_title(_('Star Chart Activity'))
                    
# Iniitialize time to now and offset to our zone.

    now = datetime.utcnow()
    (tstr, ostr) = set_time_and_UTC_offset()
    (hh, mm) = parse_zone_offset(ostr)
    zoneoffset = 60 * hh
    if (hh < 0):
      zoneoffset = zoneoffset - mm
    else:
      zoneoffset = zoneoffset + mm

# If the file StarChart.cfg exists in the Activity's data directory,
# get the longitude and latitude settings stored there.  This will
# override coordinates which might have been obtained from another
# source (e.g.: observatory.py, which is now deprecated).  (Since this
# file may not exist, the code is a bit turgid.)

    self.datafile = os.path.join(activity.get_activity_root(),\
                                 'data', 'StarChart.cfg')
    try:
      f = open(self.datafile, 'r')
    except:
      pass
    else:
      try:
        for data in f:
          if (data[0:8] == 'Latitude'):
            latitude = float(data[9:])
          elif(data[0:9] == 'Longitude'):
            longitude = float(data[10:])
      except:
        pass
      f.close()
      
# Build the translation from constellation name to constellation ID
# (needed in the "Locate" feature).

    for id in sorted(figures.keys()):
      (name, lines) = figures[id]
      abbrev_from_name[name] = id

# Create toolbox
      
    self.what_toolbar = gtk.Toolbar()
    self.where_toolbar = gtk.Toolbar()
    self.when_toolbar = gtk.Toolbar()
    self.locate_toolbar = gtk.Toolbar()
    self.about_toolbar = gtk.Toolbar()
 
    if _have_toolbox:
      toolbox = ToolbarBox()
      activity_button = ActivityToolbarButton(self)
      toolbox.toolbar.insert(activity_button, 0)
      activity_button.show()

      what_toolbar_button = ToolbarButton(
        page=self.what_toolbar,
        icon_name='toolbar-view')
      self.what_toolbar.show()
      toolbox.toolbar.insert(what_toolbar_button, -1)
      what_toolbar_button.show()

      where_toolbar_button = ToolbarButton(
        page=self.where_toolbar,
        icon_name='where')
      self.where_toolbar.show()
      toolbox.toolbar.insert(where_toolbar_button, -1)
      where_toolbar_button.show()

      when_toolbar_button = ToolbarButton(
        page=self.when_toolbar,
        icon_name='when')
      self.when_toolbar.show()
      toolbox.toolbar.insert(when_toolbar_button, -1)
      when_toolbar_button.show()

      locate_toolbar_button = ToolbarButton(
        page=self.locate_toolbar,
        icon_name='locate')
      self.locate_toolbar.show()
      toolbox.toolbar.insert(locate_toolbar_button, -1)
      locate_toolbar_button.show()

      about_toolbar_button = ToolbarButton(
        page=self.about_toolbar,
        icon_name='about')
      self.about_toolbar.show()
      toolbox.toolbar.insert(about_toolbar_button, -1)
      about_toolbar_button.show()

      separator = gtk.SeparatorToolItem()
      separator.props.draw = False
      separator.set_expand(True)
      toolbox.toolbar.insert(separator, -1)

      stop_button = StopButton(self)
      stop_button.props.accelerator = '<Ctrl>q'
      toolbox.toolbar.insert(stop_button, -1)
      stop_button.show()

      self.set_toolbar_box(toolbox)
      toolbox.show()
      self.toolbar = toolbox.toolbar
    else:
      toolbox = activity.ActivityToolbox(self)
      self.set_toolbox(toolbox)

# Fill the toolbox bars

    self._toolbar_add(self.what_toolbar, fullscreen)
    separator = gtk.SeparatorToolItem()
    separator.props.draw = True
    separator.set_expand(False)
    self._toolbar_add(self.what_toolbar, separator)
    self._toolbar_add(self.what_toolbar, button1)
    self._toolbar_add(self.what_toolbar, button2)
    self._toolbar_add(self.what_toolbar, button3)
    self._toolbar_add(self.what_toolbar, button4)
    separator = gtk.SeparatorToolItem()
    separator.props.draw = True
    separator.set_expand(False)
    self._toolbar_add(self.what_toolbar, separator)
    self._toolbar_add(self.what_toolbar, label6)
    container2.attach(rb7, 0, 1, 0, 1)
    container2.attach(rb8, 1, 2, 0, 1)
    container2.attach(rb9, 2, 3, 0, 1)
    container2.attach(rb10, 3, 4, 0, 1)
    container2.attach(rb11, 4, 5, 0, 1)
    container2.attach(rb12, 5, 6, 0, 1)
    rb7.show()
    rb9.show()
    rb11.show()
    rb8.show()
    rb10.show()
    rb12.show()
    self._toolbar_add(self.what_toolbar, container2)
    self._toolbar_add(self.where_toolbar, label1)
    self._toolbar_add(self.where_toolbar, entry1)
    container3.add(rb1)
    rb1.show()
    container3.add(rb2)
    rb2.show()
    self._toolbar_add(self.where_toolbar, container3)
    separator = gtk.SeparatorToolItem()
    separator.props.draw = False
    separator.set_expand(False)
    self._toolbar_add(self.where_toolbar, separator)
    self._toolbar_add(self.where_toolbar, label2)
    self._toolbar_add(self.where_toolbar, entry2)
    container4.add(rb3)
    rb3.show()
    container4.add(rb4)
    rb4.show()
    self._toolbar_add(self.where_toolbar, container4)
    separator = gtk.SeparatorToolItem()
    separator.props.draw = False
    separator.set_expand(False)
    self._toolbar_add(self.where_toolbar, separator)
    self._toolbar_add(self.where_toolbar, button5)
    separator = gtk.SeparatorToolItem()
    separator.props.draw = False
    separator.set_expand(False)
    self._toolbar_add(self.where_toolbar, separator)
    self._toolbar_add(self.where_toolbar, button51)
    self._toolbar_add(self.when_toolbar, rb5)
    self._toolbar_add(self.when_toolbar, rb6)
    self._toolbar_add(self.when_toolbar, label4)
    self._toolbar_add(self.when_toolbar, entry3)
    self._toolbar_add(self.when_toolbar, label5)
    self._toolbar_add(self.when_toolbar, entry4)
    separator = gtk.SeparatorToolItem()
    separator.props.draw = False
    separator.set_expand(False)
    self._toolbar_add(self.when_toolbar, separator)
    self._toolbar_add(self.when_toolbar, button6)
#    objtypecb.append_text(_('Constellations'))
    objtypecb.append_text(_('Planets'))
    objtypecb.append_text(_('Stars by Constellation'))
    objtypecb.append_text(_('Brightest Stars'))
    objtypecb.append_text(_('Deep-sky Objects'))
    self._toolbar_add(self.locate_toolbar, labell1)
    self._toolbar_add(self.locate_toolbar, objtypecb)
    (name, wbar, e, a, I, O, L0, dL) = sun
    planetscb.append_text(name)
    (name, L0, P0, N0, I, e, a, phi0, tau) = moon
    planetscb.append_text(name)
    for i in range(len(planets)):
      if (i == 2):
        pass
      else:
        (name, wbar, e, a, I, O, L0, dL) = planets[i]
        planetscb.append_text(name)
    names = []
    for code, (name, lines) in figures.iteritems():
# lines is an array of coordinates.  we ignore it.
      names = names + [name]
    for name in sorted(names):
      constscb.append_text(name)
    for i in range(len(dso_chart)):
      (nM, strCon, ra, dec, mag, majA, minA, posA, strT, strN) = dso_chart[i]
      if (strN == ''):
        dsoscb.append_text(nM)
      else:
        dsoscb.append_text(strN + ' (' + nM +')')
    container0.add(constscb)
    self._toolbar_add(self.locate_toolbar, container0)
    # container1.add(labela1)
    # labela1.show()
    # container1.add(labela2)
    # labela2.show()
    container1.add(labela3)
    labela3.show()
    container1.add(labela4)
    labela4.show()
    self._toolbar_add(self.about_toolbar, container1)
    if not _have_toolbox:
      toolbox.add_toolbar(_('What'), self.what_toolbar)
      toolbox.add_toolbar(_('Where'), self.where_toolbar)
      toolbox.add_toolbar(_('When'), self.when_toolbar)
      toolbox.add_toolbar(_('Locate'), self.locate_toolbar)
      toolbox.add_toolbar(_('About'), self.about_toolbar)

# Create the GUI objects.

    scrolled = gtk.ScrolledWindow()
    scrolled.set_policy(gtk.POLICY_NEVER, gtk.POLICY_AUTOMATIC)
    scrolled.props.shadow_type = gtk.SHADOW_NONE
    self.chart = ChartDisplay(self)
    eb = gtk.EventBox()
    vbox = gtk.VBox(False)
    self.identifyobject = gtk.Label('')
    vbox.pack_start(self.identifyobject, expand=False)
    vbox.pack_start(self.chart)
    eb.modify_bg(gtk.STATE_NORMAL, gtk.gdk.color_parse('gray'))

# Stack the GUI objects.

    scrolled.add_with_viewport(vbox)

# Connect the control widget events.

    fullscreen.connect('clicked', self._fullscreen_callback)
    button1.connect('clicked', self.chart.callback, 'night vision')
    button2.connect('clicked', self.chart.callback, 'invert display')
    button3.connect('clicked', self.chart.callback, 'flip horizontally')
    button4.connect('clicked', self.chart.callback, 'draw constellations')
    rb7.connect('clicked', self.chart.callback, 'rb7 clicked')
    rb8.connect('clicked', self.chart.callback, 'rb8 clicked')
    rb9.connect('clicked', self.chart.callback, 'rb9 clicked')
    rb10.connect('clicked', self.chart.callback, 'rb10 clicked')
    rb11.connect('clicked', self.chart.callback, 'rb11 clicked')
    rb12.connect('clicked', self.chart.callback, 'rb12 clicked')
    button5.connect('clicked', self.chart.callback, 'location change')
    button51.connect('clicked', self.chart.callback, 'home location set')
    button6.connect('clicked', self.chart.callback, 'time change')
    rb6.connect('clicked', self.chart.callback, 'user time')
    rb5.connect('clicked', self.chart.callback, 'now time')
    self.chart.connect('expose_event', self.chart.area_expose_cb)
    objtypecb.connect('changed', self.chart.callback, 'objtype sel')
    constscb.connect('changed', self.chart.callback, 'constellation sel')
    starscb.connect('changed', self.chart.callback, 'star sel')
    planetscb.connect('changed', self.chart.callback, 'planet sel')
    dsoscb.connect('changed', self.chart.callback, 'dso sel')

# Set the canvas

    self.set_canvas(scrolled)

# Show the GUI stack.

    toolbox.show()
    self.chart.show()
    eb.show()
    scrolled.show()
    self.show_all()

# FIXME: We can't do sharing yet, so hide the control for it.

    self.max_participants = 1
    if not _have_toolbox:
      toolbar = toolbox.get_activity_toolbar()
      toolbar.share.hide()

# Establish initial state of controls and do a plot.

    initialize_controls()
    self.chart.plotchart()

  def _toolbar_add(self, toolbar, component):
    item = gtk.ToolItem()
    item.add(component)
    toolbar.insert(item, -1)
    component.show()
    item.show()

  def _fullscreen_callback(self, button):
        ''' Hide the Sugar toolbars. '''
        self.fullscreen()

  def read_file(self, filename):
    global nightvision
    global invertdisplay
    global fliphorizontally
    global drawconstellations
    global limitingmagnitude
    global saved_lmag
    global latitude
    global longitude
    global specifytime
    global saved_specifytime
    global zoneoffset
    global now

    f = open(filename, 'r')
    nightvision = bool(int(self.metadata.get('Night_Vision', '0')))
    invertdisplay = bool(int(self.metadata.get('Invert', '0')))
    fliphorizontally = bool(int(self.metadata.get('Flip', '0')))
    drawconstellations = bool(int(self.metadata.get('Constellations', '1')))
    limitingmagnitude = float(self.metadata.get('Magnitude', '4.0'))
    if 'Lmag' in self.metadata:
      saved_lmag = float(self.metadata.get('Lmag'))
    specifytime = bool(int(self.metadata.get('Specify_Time', '0')))
    saved_specifytime = specifytime
    if (specifytime):
      ts = self.metadata.get('Time', now.strftime('%Y/%m/%d,%H:%M'))
      zs = self.metadata.get('Zone_Offset', str(zoneoffset))
      entry3.set_text(ts)
      entry4.set_text(zs)
      now = get_time_and_UTC_offset(entry3.get_text(), entry4.get_text())
      (hh, mm) = parse_zone_offset(entry4.get_text())
      zoneoffset = 60 * hh
      if (hh < 0):
        zoneoffset = zoneoffset - mm
      else:
        zoneoffset = zoneoffset + mm
    initialize_controls()
    self.chart.plotchart()

    
  def write_file(self, filename):
    f = open(filename, 'w')
    self.metadata['Night_Vision'] = str(int(nightvision))
    self.metadata['Invert'] = str(int(invertdisplay))
    self.metadata['Flip'] = str(int(fliphorizontally))
    self.metadata['Constellations'] = str(int(drawconstellations))
    self.metadata['Magnitude'] = str(limitingmagnitude)
    self.metadata['Latitude'] = str(latitude)
    self.metadata['Longitude'] = str(longitude)
    self.metadata['Specify_Time'] = str(int(specifytime))
# Unlike the other settings, it's easier to store the time and zone offset as
# they are represented in the text entry controls than to attempt to convert
# now and zoneoffset into a representation of local time and offset.
    self.metadata['Time'] = entry3.get_text()
    self.metadata['Zone_Offset'] = entry4.get_text()
    self.metadata['Lmag'] = str(saved_lmag)
    f.close()


  def update_config(self):
# Modify the values for latitude and longitude.
    data = []

# TODO: In future releases, StarChart.cfg may contain settings for
# other context variables.

    self.datafile = os.path.join(activity.get_activity_root(),\
                                 'data', 'StarChart.cfg')
    try:
      f = open(self.datafile, 'r')
    except:
      f = open(self.datafile, 'w')
      pass
    else:
      try:
        for old_data in f:
          if (old_data[0:8] == 'Latitude'):
            pass
          elif(old_data[0:9] == 'Longitude'):
            pass
          else:
            data[len(data)] = old_data
      except:
        pass
      f.close()
      f = open(self.datafile, 'w')
      for i in range(len(data)):
        f.write(data[i])
    f.write('Latitude=' + str(latitude) + '\n')
    f.write('Longitude=' + str(longitude) + '\n')
    f.close()
