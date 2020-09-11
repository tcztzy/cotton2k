Calibrating the model for additional sites
==========================================

This can be done by users with elementary knowledge of general agronomic
principles, soil science, meteorology, cotton physiology, and mathematics.
Some computers skills are also required.

Site-specific data files are stored in the ``\DATA\SITE`` subfolder of the
model software. It is strongly advised to make a backup of all the files in
this folder before attempting to make any changes in any of them

All files are saved as simple text files (ASCII). Notepad can be used to edit
them, but a more advanced text editor (with row and column numbering
information) is recommended.

To add a new site to the model, follow the following steps:

1. add a new line to the file ``SITELIST.DAT``. This line is in FORMAT
   (I4,1X,A20,15X,A20).

   Columns 1 to 4 contain a unique integer number, between 0 and 9999,
   assigned to this site. Do not use the same number more than once.

   Columns 6 to 25 contain the name of the site. It may be up to 20 characters
   long.

   Columns 41 to 60 contain the name of the file with the site specific
   parameters. This name may be up to 20 characters long (up to 16 characters
   plus "``.DAT``").
2. Create a new file and name it as specified in ``SITELIST.DAT``. The easiest
   way to do it is to copy an existing file (using SaveAs) and the modify it
   for the new site.

For a description of this data file go to `Structure and contents of a site data file`_

Structure and contents of a site data file
------------------------------------------

The name of this file will be not more than 20 characters long (16 characters
in the main part of the name, and “.DAT” extension added).

The first line contains the name of the site, starting in the first column. It
may be not more than 20 characters long. The next 20 lines will each contain a
site-specific parameter in the first 10 columns. Each parameter is a number
with a decimal point. Only 16 parameters are presently used, and the other
lines should be filled with zeroes (do not omit them). The model does not read
the columns from 11 to the end of the line – these columns may be used for
comments.

There are 21 lines read by the model. Additional lines may be used for comments.

The daily wind function parameters
----------------------------------

The wind speed at each hour is estimated from the daily wind run, using 4
site-specific parameters. SITEPAR(1) is the time (hours after sunrise) at
which wind speed begins to increase. SITEPAR(2) is the time (hours after solar
noon) at which wind speed is maximum. SITEPAR(3) is the time (hours after
sunset) at which wind speed reaches its night time low. SITEPAR(4) is a factor
for estimating night time wind speed.

Tip: If you have no information about wind speed at this site, use parameters from another site. If there are strong winds during the night, increase the value of SITEPAR(4).

The parameters used here and their ranges are as follows:

+-----------------+------------------+
|    parameter    |      range       |
+=================+==================+
| ``SITEPAR(1)``  | 0.0 to 2.0       |
+-----------------+------------------+
| ``SITEPAR(2)``  | 2.0 to 4.0       |
+-----------------+------------------+
| ``SITEPAR(3)``  | 0.0 to 2.5       |
+-----------------+------------------+
| ``SITEPAR(4)``  | 0.0025 to 0.0080 |
+-----------------+------------------+

Estimating the daily average dew point temperature
--------------------------------------------------

When the data for daily average dew point temperature (TDEW) are not available,
and they can not be derived from relative humidity and air temperature, they
can be crudely estimated. SITEPAR(5) is the value of TDEW when the maximum air
temperature is 20 C or lower, and SITEPAR(6) is the value of TDEW when the
maximum air temperature is 40 C or higher. Between these values of maximum air
temperature, TDEW is computed by linear interpolation.

Tip: increase the value of SITEPAR(5) if the relative humidity during the day
is higher, and increase the value of SITEPAR(6) if there are many summer
nights with dew (i.e., nightly dew point temperatures are higher than the
minimum air temperature).

The parameters used here and their ranges are as follows:

+-----------------+------------------+
|    parameter    |      range       |
+=================+==================+
| ``SITEPAR(5)``  | 6.264 to 20.0    |
+-----------------+------------------+
| ``SITEPAR(6)``  | 14.554 to 24.302 |
+-----------------+------------------+

Estimating the cloud-type correction
------------------------------------

This correction is used to decrease the computed long wave radiation emitted from the sky. The value of this correction, given by SITEPAR(7), depends on site and it is almost constant during the summer season.

The parameter used here and its range is as follows:

+-----------------+------------------+
|    parameter    |      range       |
+=================+==================+
| ``SITEPAR(7)``  | 40.0 to 120.0    |
+-----------------+------------------+

Estimating the time of daily maximum temperature
------------------------------------------------

This time at which the daily maximum air temperature occurs is used to compute
hourly air temperatures during the day. This time, in hours after solar noon,
is given by SITEPAR(8). It is site-dependant, and it is almost constant during
the summer season.

The parameter used here and its range is as follows:

+-----------------+------------------+
|    parameter    |      range       |
+=================+==================+
| ``SITEPAR(8)``  | 1.0 to 3.5       |
+-----------------+------------------+

Estimating the deep soil temperature
------------------------------------

The temperature of the last soil layer (lower boundary) is computed as a sinusoidal function of the date (Day of year), with site-specific parameters.

DPSOILT = SITEPAR(9) + SITEPAR(10) * SIN( 2.* PI * (DAYNUM - SITEPAR(11) ) / 365

The parameters used here and their ranges are as follows:

+-----------------+------------------+
|    parameter    |      range       |
+=================+==================+
| ``SITEPAR(9)``  | 20 to 24         |
+-----------------+------------------+
| ``SITEPAR(10)`` | 3.0 to 6.5       |
+-----------------+------------------+
| ``SITEPAR(11)`` | 180.             |
+-----------------+------------------+

Estimating the daily range of dew point temperature
---------------------------------------------------

The daily range of dew point temperature (TDRANGE) is computed as a linear
regression on maximum (TMAX) and minimum (TMIN) air temperatures:

``TDRANGE = SITEPAR(12) + SITEPAR(13) * TMAX + SITEPAR(14) * TMIN``

The parameters used here and their ranges are as follows:

+-----------------+------------------+
|    parameter    |      range       |
+=================+==================+
| ``SITEPAR(12)`` | -2.436 to 2.722  |
+-----------------+------------------+
| ``SITEPAR(13)`` | 0.125 to 0.820   |
+-----------------+------------------+
| ``SITEPAR(14)`` | -0.108 to –0.930 |
+-----------------+------------------+

Estimating the albedo of the soil surface
-----------------------------------------

The albedo of the soil surface is a function of its water content. When the soil surface is air-dry or drier the albedo is SITEPAR(15), and when it is at field capacity or wetter the albedo is SITEPAR(16). Between these values the albedo is computed by linear interpolation.

The parameters used here and their ranges are as follows

+-----------------+------------------+
|    parameter    |      range       |
+=================+==================+
| ``SITEPAR(15)`` | 0.12 to 0.30     |
+-----------------+------------------+
| ``SITEPAR(16)`` | 0.06 to 0.18     |
+-----------------+------------------+

The values of parameters for calibrated sites
---------------------------------------------

+------------+------------+------------+------------+------------+-----------------+
| Ca SJ West | Az Central | Il Coastal |  Il Galil  |  Il Avdat  |    Parameter    |
+============+============+============+============+============+=================+
|        0.0 |        1.0 |        1.0 |        2.0 |        0.0 | ``SITEPAR(1)``  |
+------------+------------+------------+------------+------------+-----------------+
|        3.0 |        3.0 |        4.0 |        2.0 |        2.0 | ``SITEPAR(2)``  |
+------------+------------+------------+------------+------------+-----------------+
|        2.0 |        0.0 |        2.0 |        2.5 |        1.0 | ``SITEPAR(3)``  |
+------------+------------+------------+------------+------------+-----------------+
|     0.0060 |     0.0080 |     0.0025 |     0.0025 |     0.0042 | ``SITEPAR(4)``  |
+------------+------------+------------+------------+------------+-----------------+
|      6.264 |        11. |        20. |      10.45 |     10.711 | ``SITEPAR(5)``  |
+------------+------------+------------+------------+------------+-----------------+
|     14.554 |        15. |        24. |      16.20 |    24.302  | ``SITEPAR(6)``  |
+------------+------------+------------+------------+------------+-----------------+
|        60. |        60. |       110. |       120. |       40.0 | ``SITEPAR(7)``  |
+------------+------------+------------+------------+------------+-----------------+
|        3.5 |        3.5 |        1.0 |        1.0 |        3.5 | ``SITEPAR(8)``  |
+------------+------------+------------+------------+------------+-----------------+
|        20. |        24. |        24. |        24. |        24. | ``SITEPAR(9)``  |
+------------+------------+------------+------------+------------+-----------------+
|        4.0 |        6.5 |        4.0 |        4.5 |        3.0 | ``SITEPAR(10)`` |
+------------+------------+------------+------------+------------+-----------------+
|       180. |       180. |       180. |       180. |       180. | ``SITEPAR(11)`` |
+------------+------------+------------+------------+------------+-----------------+
|      2.722 |      2.722 |     -2.436 |     -2.436 |     -2.436 | ``SITEPAR(12)`` |
+------------+------------+------------+------------+------------+-----------------+
|      0.125 |      0.125 |      0.820 |      0.820 |      0.820 | ``SITEPAR(13)`` |
+------------+------------+------------+------------+------------+-----------------+
|     -0.108 |     -0.108 |     -0.930 |     -0.930 |     -0.930 | ``SITEPAR(14)`` |
+------------+------------+------------+------------+------------+-----------------+
|       0.30 |       0.30 |       0.30 |       0.30 |       0.12 | ``SITEPAR(15)`` |
+------------+------------+------------+------------+------------+-----------------+
|       0.18 |       0.18 |       0.18 |       0.18 |       0.06 | ``SITEPAR(16)`` |
+------------+------------+------------+------------+------------+-----------------+
