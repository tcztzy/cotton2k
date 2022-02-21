# Author: Tang Ziya
# Date: 2022-02-21
subdir('ttutil')
fs = import('fs')

executable('swap',
  sources: [
    'boundbottom.f90',
    'boundtop.f90',
    'calcgrid.f90',
    'calcgwl.f90',
    'checkmassbal.f90',
    'convertdiscrvert.f90',
    'cropgrowth.f90',
    'divdra.f90',
    'drainage.f90',
    'fluxes.f90',
    'frozencond.f90',
    'functions.f90',
    'headcalc.f90',
    'hysteresis.f90',
    'initialize.f90',
    'integral.f90',
    'irrigation.f90',
    'macropore.f90',
    'macroporeoutput.f90',
    'macrorate.f90',
    'management_soil.f90',
    'meteoday.f90',
    'meteodt.f90',
    'oxygenstress.f90',
    'penmon.f90',
    'readmeteo.f90',
    'readswap.f90',
    'rootextraction.f90',
    'sharedexchange.f90',
    'sharedsimulation.f90',
    'snow.f90',
    'soilwater.f90',
    'solute.f90',
    'sptabulated.f90',
    'surfacewater.f90',
    'swap.f90',
    configure_file(
        input: '.template',
        output: 'swapoutput.utf8.f90',
        configuration: {
            'content': fs.read('swapoutput.f90', encoding: 'cp1252'),
        }
    ),
    'temperature.f90',
    'timecontrol.f90',
    configure_file(
        input: '.template',
        output: 'tridag.utf8.f90',
        configuration: {
            'content': fs.read('tridag.f90', encoding: 'cp1252'),
        }
    ),
    configure_file(
        input: '.template',
        output: 'variables.utf8.f90',
        configuration: {
            'content': fs.read('variables.f90', encoding: 'cp1252'),
        }
    ),
    'watstor.f90',
    'wofostnut.f90',
    'wofost_soil_amendments.f90',
    'wofost_soil_balancecheck.f90',
    'wofost_soil_cropresidues.f90',
    'wofost_soil_declarations.f90',
    'wofost_soil_interface.f90',
    'wofost_soil_orgmatn.f90',
    'wofost_soil_parameters.f90',
    'wofost_soil_rateconstants.f90',
    'wofost_soil_watern.f90',
  ],
  link_with: [ttutil]
)