Input Files
===========

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   creating_and_editing.rst

Use of the program follows standard WINDOWS conventions. It is assumed that the user is familiar with the common methods and conventions of using WINDOWS programs.

The input files which are used in COTTON2K are all stored as ASCII text files, and can be edited or inspected by any text editor. However, it is highly recommended to use only the facilities of the COTTON2K user interface to create, change, view or print these files.  Input data are stored in the files in metric units, but the user interface provides conversion facilities to and from English units.

To find more information about input files check the HelpForAdvancedUsers .

The following input files are used in COTTON2K:

1. Job Files.

   These files are necessary to activate one or more simulation runs.
   Each file contains in its 1st line the name of this job file, and in the next lines one or more profile file names.
   The file extension is ``*.JOB``  and they are stored in sub-folder ``..\JOBS``.

2. Profile Files.

   Each of these files defines a simulation run.
   They include information about the cultivar and the site, the specific input files to read, and the output required.
   The file extension is ``*.PRO`` and they are stored in sub-folder ``..\PROFILES``.               For more information:        Profile files

3. Climate Files.

   These files contain weather data on a daily basis.
   The file extension is ``*.ACT`` for actual weather, or ``*.PRD`` for future or predicted weather, and they are stored in sub-folder ``..\CLIMATE``.              For more information:       Weather Data files


4. Soil Hydrology Files.

   These files describe soil characteristics.
   Their file extension is ``*.HYD`` and they are stored in sub-folder ``..\SOILS``.
   For more information: Soil Hydrology files

5. Soil Initial Conditions Files.

   These files contain information on the water, nitrate and ammonium N, and organic matter in the soil at the start of the simulation.
   Their file extension is ``*.INT`` and they are stored in sub-folder ``..\SOILS``.
   For more information: Initial Soil Data Files

6. Agricultural Inputs Files.

   These files contain farm management information on irrigation, nitrogen fertilization, cultivation, defoliation, and Pix application, water table depth and salinity data, as well as data needed for predicting optimal irrigation and defoliation dates.
   Their file extension is ``*.AGI`` and they are stored in sub-folder ``..\AGINPUT``.
   For more information: Agricultural Input files

7. Plant Map Adjustment Inputs Files.

   These files contain plant mapping information, if available for the simulated fields.
   Their file extension is ``*.MAP`` and they are stored in sub-folder ``..\PLANTMAP``.
   For more information: Plant Map Adjustment files

8. Model Data Files.

   These files contain calibration data for the model.
   Their file extension is ``*.DAT`` and they are stored in sub-folders ``..\DATA``, ``..\DATA\SITE``, ``..\DATA\VARS``.
   These files should only be changed by advanced users.
   See the HelpForAdvancedUsers files for instructions.

The user interface provides the means for easily and safely creating, modifying, or deleting the job files, profile files, climate files, soil hydrology files, soil initial conditions files, agricultural inputs files, and plant map adjustments inputs files.
