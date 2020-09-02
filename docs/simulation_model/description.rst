Description of the COTTON2K Model
=================================

The main characteristics of COTTON2K, with emphasize on the differences between it and GOSSYM, will be summarized here. Most of these modifications, which have made CALGOS and COTTON2K more suitable than GOSSYM or other previous cotton models for use in the irrigated arid regions, may be summarized as follows:

1. Water Relationships.

   Potential evapotranspiration is computed on an hourly basis, using equations derived from those adopted by CIMIS (California Irrigation Management Information Service).  In order to enhance the accuracy of the computation of potential evapotranspiration, a procedure for estimating hourly values of weather parameters from the daily values has been implemented. Note that in addition to the daily weather input parameters used by GOSSYM (global radiation, maximum and minimum temperatures, rainfall, wind) daily average dew-point temperatures are now also required as input.

   The root submodel has been modified, especially concerning the response of root growth and activity to differential soil moisture conditions.  Average soil water potential is computed as an average of the root zone, weighted by root activity in each soil cell.  This soil water potential is used for computing the actual transpiration by the plants.

   Water movement in the soil is computed as a combination of implicit and explicit numerical procedures. This is done at an hourly time step.

   Leaf water potential is computed on the basis of the average soil water potential, plant resistance to water transport, and potential transpiration. The leaf water potential is then used to compute empirical water-stress factors.  These water-stress factors affect the growth rates of plant parts, aging rates of leaves and bolls, photosynthesis, and abscission rates of leaves, squares and bolls.

2. Irrigation.

   In addition to other methods of irrigation (sprinklers, furrows), the option of drip irrigation has been implemented as input to the model.

   The model can also be used to predict the irrigation requirements of the crop, under assumed weather scenarios and soil conditions, for drip as well as for other methods of irrigation.


3. Nitrogen Relationships.

   The processes of nitrogen mineralization (from decomposing organic matter) and nitrification in the soil have been modified.  Modules for denitrification, N immobilization under high C/N ratios, urea hydrolysis, and transport of nitrate and urea in the soil  have been added.

   Uptake of N by the plants is assumed to be affected by the growth requirements of each plant part, and it is simulated as a Michaelis-Menten procedure. New procedures have been devised to simulate the allocation and reallocation of nitrogen to all plant parts. Nitrogen stress factors are computed, and their effects on plant growth rates, aging of leaves and bolls, and abscission of squares and bolls are simulated.

4. Plant Growth and Phenology.

   Leaf growth is simulated separately for blades and petioles, using the monomolecular growth function. The parameters of the growth function are dependant on the node position of each leaf. The routines for leaf aging and abscission have also been completely revised.

   Boll growth is simulated separately for seed-cotton and for burrs, using improved growth functions. The logistic function is used to simulate the growth of seedcotton, whereas burr growth is assumed to be linear for the first three weeks after flowering. The routines for boll aging, abscission, and opening have also been revised.


5. Abscission of Squares and Bolls.

   The rate of abscission of squares and bolls is assumed to be affected by carbon stress, water stress, and nitrogen stress. There is a time lag (usually 5 to 6 days, depending on temperature) between the occurrence of the physiological stress and the actual abscission caused by it. The susceptibility of each square or boll to shedding is simulated as a function of its physiological age and the severity of stress.


6. Soil and Canopy Temperatures.

   The temperature of the soil surface is computed by solving the energy balance equation at the soil surface: heat conductance in the soil, incoming short wave radiation, incoming long wave radiation (from sky and from canopy), outgoing long wave radiation, sensible heat transfer, and latent heat of evaporation.

   The temperature of the plant canopy is similarly computed by solving the energy balance equation at the canopy interface with the air:  incoming short wave radiation, incoming long wave radiation (from sky and from soil), outgoing long wave radiation, sensible heat transfer, and latent heat of transpiration.

   Heat flux in the soil is computed as a combination of implicit and explicit numerical procedures. All these procedures are done at an hourly time step.

7. Time Steps used in the Model.

   Most of the procedures in the model, as in many other models, are computed in a daily time step. However, in order to increase the accuracy of the simulation, we can now utilize the enormous computing power of today's personal computers. Although weather input data are in daily values, they are used to estimate the hourly values of these data. The "heat units" (or "physiological age") used to express the effects of temperature on growth rates and phenology are now computed at an hourly time step.

   Also, the following procedures are now computed at an hourly time step: transpiration and evaporation from the soil surface, water and nitrogen movement in the soil, heat flux in the soil, energy exchanges at the soil-plant-air interfaces, soil and canopy temperatures, prediction of plant germination and emergence.

8. Use of Field Data to Modify and Adjust Model Performance.

   In some "real life" instances, the actual performance of the plants in the field may be different from the prediction of the model. This may be caused by several reasons, for example: (1) the effects of some environmental factors (like insect pests, diseases, weeds, deficiecies in P or K or micro-elements in the soil) are not simulated by the model; (2) the parameters of a new cultivar have not yet been correctly calibrated; (3) the climate and soil input variables may be unreliable.

   Many cotton growers Monitor regularly the growth and development of the plants in their fields, using a procedure called "Plant Mapping". Plant height is measured at regular intervals, as well as the number and location of fruiting sites (as squares, flowers, green bolls, mature bolls, and abscised sites). The model now allows the input of such data, computes the averages for each observation date, and creates "plant adjustment files". The model can now be rerun, using these files as input, so that its results will fit smoothly the observation data.
