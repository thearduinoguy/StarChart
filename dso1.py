# -*- coding: utf-8 -*-
# Messier Catalog
#        number.
# FIXME: Need to convert this into a dictionary such that we can more efficiently obtain
#        the astronomical data (including common name) when looking up the catalog
#        number.

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
# catalog, as they are not going  to be visible even in good binoculars.
#

from gettext import gettext as _

# Messier No.  Constellation  Right Ascension  Declination   Magnitude    Maj Axis      Min Axis    PosAngle  Type   Common Name
#     nM          strCon          dRA             dDec         dMag         dMajA         dMinA    strTyp    strName
data = [                                                  \
  ('M002',        'Aqr',       21.558333,      -0.816667,     7.50,      12.900000,    12.900000,    0.0, 'GCl', ''),  \
  ('M003',        'CVn',       13.703333,      28.383333,     7.00,      16.200000,    16.200000,    0.0, 'GCl', ''),  \
  ('M004',        'Sco',       16.393333,     -26.533333,     7.50,      26.300000,    26.300000,    0.0, 'GCl', ''),  \
  ('M005',        'Ser',       15.310000,       2.083333,     7.00,      17.400000,    17.400000,    0.0, 'GCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Butterfly_Cluster
  ('M006',        'Sco',       17.668333,     -32.216667,     4.50,      15.000000,    15.000000,    0.0, 'OCl', _('Butterfly Cluster')),  \
  # TRANS: http://en.wikipedia.org/wiki/Ptolemy's_Cluster
  ('M007',        'Sco',       17.898333,     -34.816667,     3.50,      80.000000,    80.000000,    0.0, 'OCl', _("Ptolemy's Cluster")),  \
  # TRANS: http://en.wikipedia.org/wiki/Lagoon_Nebula
  ('M008',        'Sgr',       18.063333,     -24.383333,     5.00,      60.000000,    35.000000,    0.0, 'C/N', _('Lagoon Nebula')),  \
  ('M010',        'Oph',       16.951667,      -4.100000,     7.50,      15.100000,    15.100000,    0.0, 'GCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Wild_Duck_Cluster
  ('M011',        'Sct',       18.851667,      -6.266667,     7.00,      14.000000,    14.000000,    0.0, 'OCl', _('Wild Duck Cluster')),  \
  ('M012',        'Oph',       16.786667,      -1.950000,     8.00,      14.500000,    14.500000,    0.0, 'GCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Hercules_Cluster
  ('M013',        'Her',       16.695000,      36.466667,     7.00,      16.600000,    16.600000,    0.0, 'GCl', _('Hercules Cluster')),  \
  ('M015',        'Peg',       21.500000,      12.166667,     7.50,      12.300000,    12.300000,    0.0, 'GCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Omega_Nebula
  ('M017',        'Sgr',       18.346667,     -16.183333,     7.00,      11.000000,    11.000000,    0.0, 'C/N', _('Omega Nebula')),  \
  ('M019',        'Oph',       17.043333,     -26.266667,     8.50,      13.500000,    13.500000,    0.0, 'GCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Trifid_Nebula
  ('M020',        'Sgr',       18.043333,     -23.033333,     5.00,      28.000000,    28.000000,    0.0, 'C/N', _('Trifid Nebula')),  \
  ('M021',        'Sgr',       18.076667,     -22.500000,     7.00,      13.000000,    13.000000,    0.0, 'OCl', ''),  \
  ('M022',        'Sgr',       18.606667,     -23.900000,     6.50,      24.000000,    24.000000,    0.0, 'GCl', ''),  \
  ('M023',        'Sgr',       17.946667,     -19.016667,     6.00,      27.000000,    27.000000,    0.0, 'OCl', ''),  \
  ('M025',        'Sgr',       18.526667,     -19.250000,     4.90,      40.000000,    40.000000,    0.0, 'OCl', ''),  \
  ('M028',        'Sgr',       18.408333,     -24.866667,     8.50,      11.200000,    11.200000,    0.0, 'GCl', ''),  \
  ('M030',        'Cap',       21.673333,     -23.183333,     8.50,      11.000000,    11.000000,    0.0, 'GCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Andromeda_Galaxy
  ('M031',        'And',        0.711667,      41.266667,     4.50,     178.000000,    53.000000,   34.0, 'Gal', _('Andromeda Galaxy')),  \
  # TRANS: http://en.wikipedia.org/wiki/Triangulum_Galaxy
  ('M033',        'Tri',        1.565000,      30.650000,     7.00,      73.000000,    45.000000,   22.0, 'Gal', _('Triangulum Galaxy')),  \
  ('M034',        'Per',        2.700000,      42.783333,     6.00,      35.000000,    35.000000,    0.0, 'OCl', ''),  \
  ('M035',        'Gem',        6.148333,      24.333333,     5.50,      28.000000,    28.000000,    0.0, 'OCl', ''),  \
  ('M036',        'Aur',        5.601667,      34.133333,     6.50,      12.000000,    12.000000,    0.0, 'OCl', ''),  \
  ('M037',        'Aur',        5.873333,      32.550000,     6.00,      24.000000,    24.000000,    0.0, 'OCl', ''),  \
  ('M038',        'Aur',        5.478333,      35.833333,     7.00,      21.000000,    21.000000,    0.0, 'OCl', ''),  \
  ('M039',        'Cyg',       21.536667,      48.433333,     5.50,      32.000000,    32.000000,    0.0, 'OCl', ''),  \
  ('M041',        'CMa',        6.783333,     -20.733333,     5.00,      38.000000,    38.000000,    0.0, 'OCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Orion_Nebula
  ('M042',        'Ori',        5.590000,      -5.450000,     5.00,      85.000000,    60.000000,    0.0, 'DfN', _('Orion Nebula')),  \
  # TRANS: http://en.wikipedia.org/wiki/De_Mairan's_Nebula
  ('M043',        'Ori',        5.593333,      -5.266667,     7.00,      20.000000,    15.000000,    0.0, 'DfN', _("de Mairan's Nebula")),  \
  # TRANS: http://en.wikipedia.org/wiki/Beehive_Cluster
  ('M044',        'Cnc',        8.668333,      19.983333,     4.00,      95.000000,    95.000000,    0.0, 'OCl', _('Beehive Cluster')),  \
  # TRANS: http://en.wikipedia.org/wiki/Pleiades
  ('M045',        'Tau',        3.783333,      24.116667,     1.40,     110.000000,   110.000000,    0.0, 'OCl', _('Pleiades')),  \
  ('M046',        'Pup',        7.696667,     -14.816667,     6.50,      27.000000,    27.000000,    0.0, 'OCl', ''),  \
  ('M047',        'Pup',        7.610000,     -14.500000,     4.50,      30.000000,    30.000000,    0.0, 'OCl', ''),  \
  ('M048',        'Hya',        8.230000,      -5.800000,     5.50,      54.000000,    54.000000,    0.0, 'OCl', ''),  \
  ('M050',        'Mon',        7.053333,      -8.333333,     7.00,      16.000000,    16.000000,    0.0, 'OCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Whirlpool_Galaxy
  ('M051',        'CVn',       13.498333,      47.200000,     8.00,      11.000000,     7.000000,    0.0, 'Gal', _('Whirlpool Galaxy')),  \
  ('M052',        'Cas',       23.403333,      61.583333,     8.00,      13.000000,    13.000000,    0.0, 'OCl', ''),  \
  ('M053',        'Com',       13.215000,      18.166667,     8.50,      12.600000,    12.600000,    0.0, 'GCl', ''),  \
  ('M055',        'Sgr',       19.666667,     -30.966667,     7.00,      19.000000,    19.000000,    0.0, 'GCl', ''),  \
  ('M062',        'Oph',       17.020000,     -30.116667,     8.00,      14.100000,    14.100000,    0.0, 'GCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Sunflower_Galaxy
  ('M063',        'CVn',       13.263333,      42.033333,     8.50,      10.000000,     6.000000,  104.0, 'Gal', _('Sunflower Galaxy')),  \
  ('M067',        'Cnc',        8.840000,      11.816667,     7.50,      30.000000,    30.000000,    0.0, 'OCl', ''),  \
  ('M068',        'Hya',       12.658333,     -26.750000,     9.00,      12.000000,    12.000000,    0.0, 'GCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Bodes_Galaxy
  ('M081',        'UMa',        9.926667,      69.066667,     8.50,      21.000000,    10.000000,  156.0, 'Gal', _('Bodes Galaxy')),  \
  # TRANS: http://en.wikipedia.org/wiki/Southern_Pinwheel_Galaxy
  ('M083',        'Hya',       13.616667,     -29.866667,     8.50,      11.000000,    10.000000,    0.0, 'Gal', _('Southern Pinwheel Galaxy')),  \
  ('M092',        'Her',       17.285000,      43.133333,     7.50,      11.200000,    11.200000,    0.0, 'GCl', ''),  \
  ('M093',        'Pup',        7.743333,     -23.866667,     6.50,      22.000000,    22.000000,    0.0, 'OCl', ''),  \
  # TRANS: http://en.wikipedia.org/wiki/Pinwheel_Galaxy
  ('M101',        'UMa',       14.055000,      54.350000,     8.50,      22.000000,    22.000000,    0.0, 'Gal', _('Pinwheel Galaxy')),  \
]
