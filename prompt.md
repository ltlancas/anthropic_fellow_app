I'd like you to develop a semi-analytic model for the evolution of adiabatic wind-blown bubbles as a function of radius and time. Specifically, this model should assume a spherically symmetric system with an injection of mass and energy within a set radius, and where the background medium has constant density, and solve for the fluid variables (density, pressure, temperature, and velocity) as a function of radius and time.

Using this model, along with a proper unit conversion system, I'd like you to solve for the evolution of 6 different bubbles, with input parameters given below. For each bubble, I'd like you to output the density, pressure, temperature, and velocity as a function of radius (up to 1.5* the outer radius of the bubble) at the time specified in the table. For each bubble, the energy and mass should be assumed to be injected uniformly within a sphere of radius 100 AU. All gas in each simulation should assume a mean molecular weight per hydrogen nucleus of 1.4*m_p (where m_p is the proton mass).

| Experiment Name | $\dot{M}_w$ [$M_{\odot}/{\rm yr}$] | $L_w$ [${\rm erg}\,{\rm s}^{-1}$]  | $n_{\rm H}$ [${\rm cm}^{-3}$] | $t_{\rm out}$ [${\rm kyr}$] |
| --- | --- | --- | --- | --- |
| F100 | $10^{-4}$ | $10^{38}$ | $100$ | $100$ |
| F10 | $10^{-4}$ | $10^{38}$ | $10$ | $100$ |
| F1 | $10^{-4}$ | $10^{38}$ | $1$ | $100$ |
| S100 | $10^{-2}$ | $10^{38}$ | $100$ | $100$ |
| S10 | $10^{-2}$ | $10^{38}$ | $10$ | $100$ |
| S1 | $10^{-2}$ | $10^{38}$ | $1$ | $100$ |

For each of the experiments there should be two output files, the first (`<experiment_name>.csv`) should specify the fluid variables over the radius range specified above. The second (`<experiment_name>.txt`) should specify two numbers, the radius of the bubble in parsecs and the pressure of the hot gas within the bubble (divided by Boltzmann's constant) in units of Kelvin per cubic centimeter.

You may reference the scientific literature in order to help develop this model. In developing the model I'd like you to state any assumptions that you've made in an `assumptions.md` file. This file should include any derivations that you carried out in building the model. The semi-analytic model should be built as a python library and can make use of numpy and/or scipy (make sure these libraries are installed locally but in a standalone environment meant only to be used for the purposes of this task).

Once this is completed I would like you to develop a one-dimensional, finite-volume hydrodynamics code that works in spherical polar coordinates. Implement versions that work at first and second order in both time and space. Record documentation on the structure of the code and how to run it in a `code_documentation.md` file.

Then run this code on the above F100 problem and use the semi-analytic model you developed to inform the size of the domain. Run the simulation at first and second order and 4 different spatial resolutions where each increases the number of resolving elements by a factor of 2 from the last. Write the outputs in a `F100_sim_VV_NXXXXX_snapYYY.csv` where `VV` is either `fo` or `so` for first or second order, `XXXXX` specifies the number of resolution elements in the simulation, `YYY` specifies the snapshot number, and the file contains all primitive fluid variables and the radial position of each grid point. There should also be a short `F100_sim_NXXXXX_times.txt` text file that contains the times of each output. There should be 5 different snapshots, each equally spaced in time.

For each simulation snapshot, use your semi-analytic model to calculate each of the fluid variables at each of the grid points and perform an L2 norm comparison between the semi-analytic model and the simulation output. Store the average L2 norm for each fluid variable (density, temperature, pressure, and velocity) in separate `ZZZZZ_L2_norm.txt` where `ZZZZZ` is one of the fluid variables and each such file has one column containing the resolution, one column containing the average L2 norm for the first order, and the last column containing the average L2 norm for the second order simulations.
