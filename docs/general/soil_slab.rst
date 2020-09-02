The Simulated Soil Slab
=======================

All soil processes are simulated in an average soil slab placed perpendicular to the direction of the plant rows. The thickness of this slab is 1 cm. Its width is the same as the average spacing between plant rows, and its depth is 200 cm from the soil surface.

When there is a uniform distance between plant rows, it is assumed that the plant is situated in the middle of the slab. All distances are measured from the left edge of the slab. Thus, if the distance between rows is 96 cm, the plant is 48 cm from the edge of the slab.

The neighboring slab is like a mirror image of the defined slab. Thus, if the horizontal distance of a drip irrigation line is input as 43 cm, it means that it is situated 5 cm from each plant row. If its distance is input as 0 (zero) cm, it means that there are drip lines for each two rows of plants, placed in the middle of the space between these rows.

When "Skip-rows" are used, the user defines skip row width as the smaller distance between two adjacent rows.  In this case, "row spacing" in the input will mean the larger distance between rows. The slab width will be the average of these numbers. Note that in all model output files, and in all model computations, "row spacing" will mean the average spacing (which is the same as slab width).


.. figure:: /_static/slabreg.png
   :width: 256px
   :height: 320px
   :alt: Slab regular

   Fig. 1 Regular spacing, a drip line for each plant row.

.. figure:: /_static/slabskip.png
   :width: 256px
   :height: 320px
   :alt: Slab skip

   Fig. 2 Skip rows, one drip line for two plant rows.
