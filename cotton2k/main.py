"""
The main user interface of the Cotton2K model.
"""
import PySimpleGUI as sg

# sg.theme("darkAmber")


# All the stuff inside your window.
layout = [
    [
        sg.Menu(
            [
                ["&File", ["&Open Job Files::-OPEN-JOB-", "E&xit"]],
                [
                    "Edit Input Files",
                    [
                        "Profiles",
                        "Initial Soil Data",
                        "Soil Hydrology",
                        "Agricultural Inputs",
                        "Weather Data",
                        "Plant Map Adjustment",
                    ],
                ],
                ["Check Output Files", ["Files", "Soil Maps", "Charts"]],
                ["&Help", ["&Help Topics", "---", "&About Cotton2K"]],
            ]
        )
    ],
    [sg.Text("Job File Name:"), sg.Input(disabled=True)],
    [
        sg.Column(
            [
                [sg.Text("AVAILABLE PROFILES:")],
                [sg.Listbox(values=[], size=(30, 40))],
            ],
        ),
        sg.Column(
            [
                [sg.Button("Add Profile", key="-ADD-PROFILE-")],
                [sg.Button("Remove Profile")],
                [sg.Button("Save Job")],
                [sg.Button("Run")],
            ]
        ),
        sg.Column(
            [
                [sg.Text("RUN QUEUE:")],
                [sg.Listbox(values=[], size=(30, 40))],
            ],
        ),
    ],
]

# Create the Window
window = sg.Window("Cotton2K - Cotton Simulation Model - Main User Interface", layout)
# Event Loop to process "events" and get the "values" of the inputs
while True:
    event, values = window.read()
    if event in (sg.WIN_CLOSED, "Exit"):  # if user closes window or clicks cancel
        break
    if event == "-ADD-PROFILE-":
        pass
    print("You entered ", values[0])
    print("Event:", event)

window.close()
