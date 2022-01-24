# PoissonFDM

Solve the generalized Poisson equation using the finite-difference method (FDM), following [Dr. James Nagel's lecture](https://my.ece.utah.edu/~ece6340/LECTURES/Feb1/Nagel%202012%20-%20Solving%20the%20Generalized%20Poisson%20Equation%20using%20FDM.pdf).


Initially, port James' MATLAB code for solving 2D systems:
http://www.drjamesnagel.com/codes/Poisson_FDM_Solver_2D.m

TODO:
1. Finish Julia port of MATLAB code.  Key sticking points:
  a. Matrix inversion (MATLAB: V = A\b;)
  b. Switch/case statement (Julia: @match)
  c. Sparse array construction
2. Get example charge systems (e.g. parallel-plate capacitor) working with visualization.
3. Investigate runtime vs. system size.
  --> check visualization render time as a proportion of total runtime, too!
4. Expand to 3D systems.
