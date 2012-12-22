# -*- coding: utf-8 -*-
# NGC Catalog subset
# FIXME: Need to convert this into a dictionary such that we can more efficiently obtain
#        the astronomical data (including common name) when looking up the catalog
#        number.
#
#     Key to Type representation:
#	SNR = supernova remnant
#	GCl = globular cluster
#	OCl = open cluster
#	C/N = cluster and nebula
#	PlN = planetary nebula
#       Gal = Galaxy
#	Dbl = double star
#	DfN = Diffuse Nebula
#	??? = unknown or unclassified object
#
# Note: Objects smaller than 10' or dimmer than magnitude 9.0 have been removed from this
# catalog, as they are not going  to be visible even in good binoculars.  Objects which
# duplicate entries in the Messier catalog have also been removed.
#
# Most of these objects come from the Bennett list of southern sky objects which are not
# comets.  The rest are NGC objects which have names.
#
#

from gettext import gettext as _

#   NGC No.    Constellation  Right Ascension  Declination   Magnitude    Maj Axis      Min Axis    PosAngle  Type   Common Name
#     nM          strCon          dRA             dDec         dMag         dMajA         dMinA    strTyp    strName
data = [                                                  \
#
# Need NGC number for an object close to the LMC center coordinates 
#
  ('LMC',       'Dor',        5.392913,     -69.756111,    0.90,     645.000000,   550.200000,    0.0, 'Gal', _('Large Magellanic Cloud')),  \
#
#
#
  ('NGC0055',       'Scl',        0.248333,     -39.183333,    7.87,      32.400000,     5.600000,    0.0, 'Gal', ''),  \
  ('NGC0104',       'Tuc',        0.401666,     -72.083333,    4.91,      30.900000,    30.900000,    0.0, 'GCl', _('47 Tucanae')),  \
  ('NGC0253',       'Scl',        0.793333,     -25.283333,    8.00,      27.500000,     6.800000,    0.0, 'Gal', _('Sculptor Galaxy')),  \
  ('NGC0288',       'Scl',        0.880000,     -26.583333,    8.10,      13.800000,    13.800000,    0.0, 'GCl', ''),  \
  ('NGC0292',       'Tuc',        0.879111,     -72.828611,    2.70,     320.000000,   185.000000,    0.0, 'Gal', _('Small Magellanic Cloud')),  \
  ('NGC0300',       'Scl',        0.914999,     -37.683333,    9.00,      21.900000,    15.500000,    0.0, 'Gal', ''),  \
  ('NGC0362',       'Tuc',        1.053333,     -70.849999,    6.40,      12.900000,    12.900000,    0.0, 'GCl', ''),  \
  ('NGC0884',       'Per',        2.373333,      57.116666,    4.30,      60.000000,    60.000000,    0.0, 'OCl', _('Double Cluster')),  \
  ('NGC1291',       'Eri',        3.288333,     -41.133333,    8.50,      10.500000,     8.100000,    0.0, 'Gal', ''),  \
  ('NGC1499',       'Per',        4.011666,      36.616666,    6.00,     145.000000,    30.000000,    0.0, 'DfN', _('California Nebula')),  \
#
# The LMC overlaps NGC2070, so this object is not useful in the catalog.
#
#  ('NGC2070',       'Dor',        5.643333,     -69.083333,    8.20,      40.000000,    20.000000,    0.0, 'DfN', 'Tarantula Nebula'),  \
#
  ('NGC2264',       'Mon',        6.685000,       9.883333,    3.90,      20.000000,    20.000000,    0.0, 'C/N', _('Cone Nebula')),  \
  ('NGC2267',       'Pup',        7.876666,     -26.383333,    7.00,      16.000000,    16.000000,    0.0, 'C/N', ''),  \
  ('NGC2237',       'Mon',        6.516666,       4.949999,    9.00,      78.000000,    78.000000,    0.0, 'C/N', _('Rosette Nebula')),  \
  ('NGC2627',       'Pyx',        8.621666,     -29.949999,    8.00,      11.000000,    11.000000,    0.0, 'C/N', ''),  \
  ('NGC2808',       'Car',        9.200000,     -64.866666,    6.30,      13.800000,    13.800000,    0.0, 'GCl', ''),  \
  ('NGC3201',       'Vel',       10.293333,     -46.416666,    6.80,      18.200000,    18.200000,    0.0, 'GCl', ''),  \
  ('NGC3242',       'Hya',       10.413333,     -18.633333,    9.00,      20.800000,    20.800000,    0.0, 'PlN', _('Ghost of Jupiter')),  \
  ('NGC3372',       'Car',       10.729999,     -59.866666,    1.00,     120.000000,   120.000000,    0.0, 'DfN', _('Eta Carina  Nebula')),  \
  ('NGC4372',       'Mus',       12.429999,     -72.666666,    7.80,      18.600000,    18.600000,    0.0, 'GCl', ''),  \
  ('NGC4755',       'Cru',       12.893333,     -60.333333,    4.20,      10.000000,    10.000000,    0.0, 'OCl', _('Jewel Box')),  \
  ('NGC4833',       'Mus',       12.993333,     -70.883333,    7.40,      13.500000,    13.500000,    0.0, 'GCl', ''),  \
  ('NGC4945',       'Cen',       13.090000,     -49.466666,    9.00,      20.000000,     3.800000,    0.0, 'Gal', ''),  \
  ('NGC5128',       'Cen',       13.424999,     -43.016666,    7.00,      25.700000,    20.000000,    0.0, 'Gal', _('Centaurus A')),  \
  ('NGC5139',       'Cen',       13.446666,     -47.483333,    3.70,      36.300000,    36.300000,    0.0, 'GCl', _('Omega Centauri')),  \
  ('NGC5617',       'Cen',       14.496666,     -60.716666,    6.30,      10.000000,    10.000000,    0.0, 'OCl', ''),  \
  ('NGC5897',       'Lib',       15.289999,     -21.016666,    8.60,      12.600000,    12.600000,    0.0, 'GCl', ''),  \
  ('NGC5927',       'Lup',       15.466666,     -50.666666,    8.30,      12.000000,    12.000000,    0.0, 'GCl', ''),  \
  ('NGC6362',       'Ara',       17.531666,     -67.050000,    8.30,      10.700000,    10.700000,    0.0, 'GCl', ''),  \
  ('NGC6397',       'Ara',       17.678333,     -53.666666,    5.70,      25.700000,    25.700000,    0.0, 'GCl', ''),  \
  ('NGC6541',       'CrA',       18.133333,     -43.699999,    6.60,      13.100000,    13.100000,    0.0, 'GCl', ''),  \
  ('NGC6723',       'Sgr',       18.993333,     -36.633333,    7.30,      11.000000,    11.000000,    0.0, 'GCl', ''),  \
  ('NGC6744',       'Pav',       19.163333,     -63.849999,    9.00,      20.000000,    12.900000,    0.0, 'Gal', ''),  \
  ('NGC6752',       'Pav',       19.181666,     -59.983333,    5.40,      20.400000,    20.400000,    0.0, 'GCl', ''),  \
  ('NGC6822',       'Sgr',       19.748333,     -14.800000,    9.00,      15.500000,    13.500000,    0.0, 'Gal', _("Barnard's Galaxy")),  \
  ('NGC6888',       'Cyg',       20.200000,      38.350000,    7.40,      18.000000,    12.000000,    0.0, 'DfN', _('Crescent Nebula')),  \
  ('NGC6960',       'Cyg',       20.761666,      30.716666,    7.00,     180.000000,   180.000000,    0.0, 'DfN', _('Veil Nebula')),  \
  ('NGC7000',       'Cyg',       20.979999,      44.333333,    4.00,     120.000000,   100.000000,    0.0, 'DfN', _('North America Nebula')),  \

]
