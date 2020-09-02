.. cotton2k documentation master file, created by
   sphinx-quickstart on Tue Sep  1 20:14:49 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
============

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Model description </simulation_model/description>
   /general/validation_and_calibration
   Upgrading from Previous Versions </general/upgrading_from_previous_versions.rst>
   /general/soil_slab
   /input_files/index

Welcome to the COTTON2K cotton simulation model,  version 4.0, released September 2004.

The COTTON2K cotton simulation model has been originally developed in 1992 from GOSSYM-COMAX.  Its main purpose was to make the model more useful for simulating cotton production under irrigation in the arid regions of the Western US and in Israel.  Since many changes have been made in the model, it has been given a new name: CALGOS (for CALifornia GOSsym).  CALGOS, like GOSSYM, is a process-level model.  It simulates the processes occurring  in the soil, plant, and the microenvironment, and the interactions between these processes and the management inputs applied to the field. The year 2000 version has been renamed COTTON2K.

The 1997 version of CALGOS had options both for Windows 3.1 and for Windows 9x. The user interface had been compiled by Microsoft Visual C++ version 4, and the model itself by Microsoft PowerStation Fortran version 4.

The 1998 and 1999 versions of CALGOS could be run by Windows 95, 98 and Windows NT4 only. The user interface had been compiled by Microsoft Visual C++ version 5, and the model itself by Digital Visual Fortran version 5. The COTTON2K user interface has been compiled by Microsoft Visual C++ version 6. It makes full use of Windows 9x/NT conventions for editing, opening and saving files, etc. The model itself has been compiled by Compaq Visual Fortran version 6.5, and it runs in the Windows environment.

Version 1.1 of Cotton2K fixed some bugs and problems that have been brought to our attention in the previous versions. The date format, which had been previously in the USA style (mm/dd/yy), has been changed to a more universally accepted format (dd-MON-yyyy, as in the example "12-JUN-2001"). The use of a 4-digit year has the added advantage of resolving all the turn of the century problems caused by a 2-digit year specification. This version also enables the users to create their own cultivar-specific and site-specific data files.

After the release of Version 1.1 of Cotton2K many procedures, mainly for soil and canopy temperature, soil water and soil nitrogen, climate data handling, and the response of the plants to stresses, have been rewritten and enhanced. Also were added capabilities to handle the effects of high water table and salinity.

The user interfaces for input and output have been made more user-friendly, and recompiled using the C++ compiler of Microsoft's Visual.Net 2003 system. The code of the model itself has been translated from Fortran to C++, and has also been compiled by Visual.Net 2003.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
