
THE MAIN USER INTERFACE SCREEN
==============================

The COTTON2K program can be executed from the Desktop by the standard Windows methods.  The main user interface appears on the screen when the program starts.

Before attempting to run a simulation, you have to prepare the necessary input files. The PROFILE FILE defines a simulation run and contains the names of other input files to be used.  Each simulation is uniquely defined in it.  Several simulations may be run one after the other, and their run sequence is defined in the JOB FILE , which contains the names of one or more Profile Files.  To run the COTTON2K model, a Job file must first be opened or created.

The user interface consists of the following elements:  a menu, file lists, boxes, and buttons.

The Menu
    FILE menu
        Open Job Files
            Open an existing or a new JOB file (you can also use ctrl / O for this purpose).
        Exit
            Exit program (you can also use alt / F4 for this purpose).
    EDIT INPUT FILES menu
        Profiles
            Creating or editing Profile input files.
        Initial Soil Data
            Creating or editing Soil Initial Conditions input files.
        Soil Hydrology
            Creating or editing Soil Hydrology data input files.
        Agricultural Inputs
            Creating or editing Agricultural input files.
        Weather Data
            Creating or editing Weather Data input files.
        Plant Map Adjustment
            Creating or editing Plant Map Adjustment input files.
    CHECK OUTPUT FILES menu
        Files
            Examine existing output files as texts.
        Charts
            Examine charts from output files.
        Soil Maps
            Examine soil maps from output files.
    MAKE PLANT MAP FILES
          Opens a special program to record and process plant observations in the field. See ????to add href???? for more information.
    HELP menu
        Help Topics
            Use this HELP file.
        About COTTON2K
            Information about this program.

File Lists
    Available Profiles
        List of all available profile files.

    Run Queue
        List of profile files selected to be run in this Job.

Edit Box
          Job File Name          Contains the name of the current Job File.

Buttons

    The following buttons are shown:

        <Add Profile>

        <Remove Profile>

        <Save Job>

        <Run>

    Each button may sometimes be disabled, and in this case it will usually be pale gray (or have the typical color defined by your Windows Desktop).  A button is activated by clicking on it with the left button of the mouse, but it will not respond when disabled.

    For more information:

        Using the Job File Dialog

        Running a Job

        Creating and Editing Input Files
