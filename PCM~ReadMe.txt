Zanil Narsing
PCM Model

Starting on line 4, the constants and initial conditions are defined
	The values for the material were found using source 4 (Zabla, Marin, Cabeza, Mehling, 2003).
	The values are for Paraffin Wax. This should be updated upon further study. The thermophysical properties for Parrifin are easily found.
	
Starting on line 35, the initial conditions for the temperature node is defined
	The initial temperature of the PCM is 290 K
	The temperature of the surroundings is also 290 K
	The melting point of the PCM is 315 K.

Starting on line 39, the temperature profile of the outter pipe wall is defined. This is the only source of heat in the model.
	This model includes a charge time, t_charge, which can be defined to fit with the steam acc. model.
	T_max is the maximum temperature reached during charging, and the temperature held during storage

The for loop starting on line 45 is used to set up a linear charging time for the wall temperature.

Lines 54, 55, and 56 establish initial conditions for the temperature array.
	At the first radial location (wall location), the temperature is always the temperature of the wall.
	At time (t=0), all radial locations are at T=T_initial

Lines 62-75 establish arrays and vectors for useful variables
	R_pc is the radial location of the interface between solid and liquid
	R_melt is the radial location of the interface between temperatures above the melting point and temperatures below the melting point
	The distance between R_pc and R_melt is fixed in this model (constant boundary layer thickness)

Line 77 initializes ifmelted at 0.
	ifmelted determines where the code draws previous time step temperature profiles
	previous time step profiles are needed when using the finite difference method
	if ifmelted == 0, then the code can draw temperatures from a non-melted time step
	if ifmelted == 1, then the code can draw temperatures from a melted time step

Line 81 begins the time loop
	this loop steps through time with the index "m"

Line 82 checks if the ifmelted is 0 or 1
	if 0, the PCM has not yet melted, and can procede with solid-only diffusion
	
The for loop from lines 85 through 96 are the solid-only space loop.
	This loop incorporates solid heat transfer using the finite difference method. This equation can be found in source 1, "Heat Conduction in Cylindrical and Spherical Coordinates II"

Line 97 establishes the Temp_check array.
	This array is used to determine if part of the PCM is above the melting temperature, and thus suited for the phase change code

Line 99 begins the if statement used to determine if the PCM is melting

ifmelted is turned to 1 in line 101

The if statement on starting on line 103 is used to correct for initial conditions
	if the index for the solid to liquid interface (n_pc) is exactly 0, then it is changed to equal the index of Temp_max (the greatest temperature within Temp_check)
	R_pc is then set equal to r at the correct location
	R_melt is then set to be 1 radial step after R_pc
	if the solid to liquid interface is not exactly 0, then Qdot_in, Qdot_out, and Qdot_store need to be calculated. 
	The equations for Qdot_in and Qdot_out can be found in source 2, (Incropera, 2011). (This is the Heat Transfer Text book)
	The value Qdot_store is used to move the interface location (R_pc)

Lines 130-150 find the temperature profiles in the 3 regions
	Region 1: solid region - solid heat diffusion
		from 1 radial location past R_melt to the end of the PCM
	Region 2: liquid region - liquid heat diffusion
		from the inner wall to one radial location before R_pc
	Region 3: boundary layer - Temperature held at melting temperature
		includes nodes of R_pc and R_melt

The ifmelted if statement continues on line 154 with the else statement
	This is run if the ifmelted is not equal to 0, meaning, the code already ran a case where it established an R_pc and an R_melt
	This set of equations is exactly the same as the in the previous section
	The reason for this section of code is to draw temperatures from the previous time step without re-establishing initial conditions

Lines 170-194 find the temperature profiles in the 3 regions in the same way as before

The for loops end on line 196

A plot is created on line 201, in which it shows the temperature profile from many different time steps



