In answering the tasks laid out below **do not** read or consider any information laid out in `grader.md`, `grader_tools/`, `exp_*/`, or a memory directory associated with this project. This task should take less than 15 minutes.

# Output File Formats
All output paths given below are relative to the directory you are working in. The output files will be read by an automated script, so follow these formats exactly:

- **`.csv` files** (`sam_outputs/<experiment_name>.csv` and `sim_outputs/F100_sim_VV_NXXXXX_snapYYY.csv`): comma-separated. The first line must be exactly `r,rho,v,P,T`, with no leading `#` and no other header or comment lines. Each following line gives the values at one radius, in these columns and units:

| Column name | Quantity | Units |
| --- | --- | --- |
| `r` | radius | ${\rm pc}$ |
| `rho` | mass density | ${\rm g}\,{\rm cm}^{-3}$ |
| `v` | radial velocity | ${\rm km}\,{\rm s}^{-1}$ |
| `P` | pressure divided by Boltzmann's constant ($P/k_B$) | ${\rm K}\,{\rm cm}^{-3}$ |
| `T` | temperature | ${\rm K}$ |

- **`sam_outputs/<experiment_name>.txt`**: a single line with exactly two numbers separated by a space: the outer radius of the bubble in ${\rm pc}$, then the hot-gas pressure $P/k_B$ in ${\rm K}\,{\rm cm}^{-3}$. No labels, units, or comments.
- **`sim_outputs/F100_sim_NXXXXX_times.txt`**: one number per line, where line $k$ is the time of snapshot $k$ in ${\rm kyr}$. No header or comments.
- **`convergence/<variable>_L1_norm.txt`**, where `<variable>` is exactly one of `density`, `temperature`, `pressure`, or `velocity`: one line per resolution, in order of increasing resolution. Each line has exactly three numbers separated by spaces: the number of cells $N$, the average L1 norm of the first-order simulation, and the average L1 norm of the second-order simulation. No header or comments.
- **Numbers in file names:** `VV` is `fo` (first order) or `so` (second order). `XXXXX` is the number of cells, zero-padded to 5 digits (e.g. `N01024`). `YYY` is the snapshot number, zero-padded to 3 digits and starting at `001`.
- Numbers inside files may use plain decimal or scientific notation (e.g. `7.63` or `7.63e+00`).

# Task 1: Semi-Analytic Model
I'd like you to develop a semi-analytic model for the evolution of adiabatic wind-blown bubbles as a function of radius and time. Specifically, this model should assume a spherically symmetric system with an injection of mass and energy within a set radius, and where the background medium has constant density, and solve for the fluid variables (density, pressure, temperature, and velocity) as a function of radius and time.

Using this model, along with a proper unit conversion system, I'd like you to solve for the evolution of 3 different bubbles, with input parameters given below. For each bubble, I'd like you to output the density, pressure, temperature, and velocity as a function of radius (up to 1.5* the outer radius of the bubble) at the time specified in the table. For each bubble, the energy and mass should be assumed to be injected uniformly within a sphere of radius 100 AU. All gas in each simulation should assume a mean molecular weight per hydrogen nucleus of 1.4*m_p (where m_p is the proton mass).

| Experiment Name | $\dot{M}_w$ [$M_{\odot}/{\rm yr}$] | $L_w$ [${\rm erg}\,{\rm s}^{-1}$]  | $n_{\rm H}$ [${\rm cm}^{-3}$] | $t_{\rm out}$ [${\rm kyr}$] |
| --- | --- | --- | --- | --- |
| F100 | $10^{-4}$ | $10^{38}$ | $100$ | $100$ |
| F10 | $10^{-4}$ | $10^{38}$ | $10$ | $100$ |
| F1 | $10^{-4}$ | $10^{38}$ | $1$ | $100$ |

For each of the experiments there should be two output files, the first (`sam_outputs/<experiment_name>.csv`) should specify the fluid variables over the radius range specified above. The second (`sam_outputs/<experiment_name>.txt`) should specify two numbers, the outer radius of the bubble (i.e. the location of the outer shock) in parsecs and the pressure of the hot gas within the bubble (divided by Boltzmann's constant) in units of Kelvin per cubic centimeter. Both files should follow the formats in "Output File Formats" above.

You may reference the scientific literature in order to help develop this model. In developing the model I'd like you to state any assumptions that you've made in an `sam_outputs/assumptions.md` file. This file should include any derivations that you carried out in building the model. The semi-analytic model should be built as a python library and can make use of numpy and/or scipy (make sure these libraries are installed locally but in a standalone environment meant only to be used for the purposes of this task). Once this task is completed, all code created for this task should be set and not edited for the rest of the prompt, including after the partial or whole completement of the next task.

# Task 2: One-dimensional Hydrodynamics Code
Once the semi-analytic model is developed I would like you to develop a one-dimensional, finite-volume hydrodynamics code that works in spherical polar coordinates. Implement versions that work at first and second order in both time and space. Record documentation on the structure of the code and how to run it in a `code_documentation.md` file. You may develop tests for the code before performing the tests laid out below, but once you have moved on to running the below simulations you may not further modify the code. The code tests should be designed to run in less than 10 seconds.

Then run this code on the above F100 problem (and **only** on this problem) but with an injection radius of 0.5 pc, and use the semi-analytic model you developed to inform the size of the domain. Run the simulation at first and second order and 3 different spatial resolutions where each increases the number of resolving elements by a factor of 2 from the last. Run each simulation from $t = 0$ to $t = 100\,{\rm kyr}$ and write 5 snapshots, at $t = 20, 40, 60, 80,$ and $100\,{\rm kyr}$, numbered `001` to `005`. Do not write a snapshot of the initial conditions. Write the outputs in a `sim_outputs/F100_sim_VV_NXXXXX_snapYYY.csv` where `VV` is either `fo` or `so` for first or second order, `XXXXX` specifies the number of resolution elements in the simulation, `YYY` specifies the snapshot number, and the file contains all primitive fluid variables and the radial position of each grid point. There should also be a short `sim_outputs/F100_sim_NXXXXX_times.txt` text file that contains the times of each output, one per line in order of snapshot number, in units of ${\rm kyr}$. All of these files should follow the formats in "Output File Formats" above. The simulation code should not be modified at all after you have moved on to Task 3.

# Task 3: Comparison and Convergence
For the last simulation snapshot, use your semi-analytic model to calculate each of the fluid variables at each of the grid points and perform an L1 norm comparison between the semi-analytic model and the simulation output. The L1 norm should be computed from relative errors, $|q_{\rm sim} - q_{\rm SAM}|/|q_{\rm SAM}|$ for each fluid variable $q$, averaged over all grid points with radii smaller than the outer radius of the bubble (as given by the semi-analytic model at that time). Store the average L1 norm for each fluid variable (density, temperature, pressure, and velocity) in separate `convergence/<variable>_L1_norm.txt` files, where `<variable>` is one of `density`, `temperature`, `pressure`, or `velocity`. Each file has one column containing the resolution (the number of cells $N$), one column containing the average L1 norm for the first order, and the last column containing the average L1 norm for the second order simulations, following the format in "Output File Formats" above.

Record the total time this task took in a `time.txt` file.
