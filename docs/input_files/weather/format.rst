FORMAT OF WEATHER INPUT FILES
=============================

Both Actual ( *.ACT ) and Predicted ( *.PRD ) weather file types have the same basic structure:

The first line of the file :

    columns 1 - 30       contain optional text describing the file (location, year, etc.)

    column 33              indicates the units used for radiation:  0 = Langleys;  1 = MJ / m2;  2 = W/m 2

    column 36              indicates units used for minimum and maximum temperatures:  0 = Fahrenheit;  1 = Centigrade

    column 39              indicates units used for rain:  0 = inches;  1 = mm

    columns 41 - 42      indicate units used for wind run:  0 = miles per day;  1 = km per day;  -1 = data not available (in this case, input of seasonal or estimated average is required)

    columns 44 - 45     indicate units used for dew-point temperatures:  0 =Fahrenheit;  1 = Centigrade;  -1 = data not available

    columns 61 - 70     seasonal average of wind run, km per day (this should be supplied if columns 41- 42 are -1).

    columns 77 - 80     the (4 digit) year of the first day in the file.  Note: The year will be 0 or blank for *.PRD files.

All other lines of the file :

    Each line contains daily values for each day, in the following order:

    (1) The Julian date (days from the beginning of the year)

    (2) The date string (or a blank string for *.PRD files) enclosed in single quotes.

    (3) Radiation.

    (4) Maximum temperature.

    (5) Minimum temperature.

    (6) Rainfall.

    (7) Daily wind run (leave blank if not available - when columns 41 - 42 of the first line are -1).

    (8) Dew point temperature (leave blank if not available - when columns 44 - 45 of the first line are -1).

An example of the first lines of a file :



Westside CIMIS weather 1994     0  1  1  1  1                     0.00      1994

  60  '01-MAR-1994'   440.11  24.20   6.20  0.000 141.00  10.85

  61  '02-MAR-1994'   438.02  25.35   8.55  0.000 136.50  11.20

  62  '03-MAR-1994'   423.65  25.58   9.50  0.000 194.80  11.10



Notes

Although COTTON2K can read the climate file (with the exception of the first line) in a free format, EDCLIM2K  creates and expects a fixed format for the data lines, as follows:

    (1) The Julian date (day of year)  - in columns 2 - 4.

    (2) The date string (including the single quotes) - in columns 7 - 19.  The file will have the date as 'dd-MON-yyyy'.

    (3) Radiation - in columns 22 - 28.

    (4) All the other data - in 5 successive 7-column blocks (in columns 29 to 63).

    All dates in the simulation period must have weather data, from at least one weather file (actual or predicted).

Although COTTON2K can read the data in English or metric units, all data are saved by EDCLIM2K  in metric units. However, the user may optionally use the data in English units in the spreadsheet, charts or printed documents.