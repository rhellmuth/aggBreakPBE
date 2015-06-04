Aggregation-Breakup Equation
============================

The Python code **aggBreakPBE.py** is used to simulate the aggregation and breakup of colloid particles using the Smoluchowski coagulation equation for aggregation and the empirical breakup model of Pandya and Spielman. This code was originally designed for simulating the aggregation of blood platelets, but it can be generalized for any kind of colloid particle without much trouble. Be aware that the code may contain bugs. Please don't hesitate in contacting me if you find any bug.

Aggregation rate is related with the statistical average collision of particles. The collision can occur by Brownian motion, shear flow (gradient of velocity), or by an enhanced diffusion effect caused by the tumbling and rolling of red blood cells under motion. The Sorensen enhanced diffusivity model is used for that.

The class **aggBreak** possesses methods for setting the system up, solving it, processing the data, and saving the simulation case in cvs files for further analysis.

References
==========

1. von Smoluchowski, M. R. Versuch einer mathematischen Theorie der Koagulationskinetik kolloider Lösungen. Zeitschrift fuer Phys. Chemie 92, 129–168 (1917).

2. Pandya, J. D. & Spielman, L. A. Floc breakage in agitated suspensions: Effect of agitation rate. Chem. Eng. Sci. 38, 1983–1992 (1983).

3. Kumar, S. & Ramkrishna, D. On the solution of population balance equations by discretization—I. A fixed pivot technique. Chem. Eng. Sci. 51, 1311–1332 (1996).

4. Garrick, S. C., Lehtinen, K. E. J. & Zachariah, M. R. Nanoparticle coagulation via a Navier–Stokes/nodal methodology: Evolution of the particle field. J. Aerosol Sci. 37, 555–576 (2006).

5. Bäbler, M. U. & Morbidelli, M. Analysis of the aggregation-fragmentation population balance equation with application to coagulation. J. Colloid Interface Sci. 316, 428–41 (2007).
