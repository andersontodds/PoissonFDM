## PoissonFDM.jl
# Todd Anderson
# Jan 24 2022
#
# Solve the generalized Poisson Equation using the finite-difference method (FDM).
# Adapted from https://drjamesnagel.com/codes/Poisson_FDM_Solver_2D.m
# Further reading:
#   Nagel, J. R, "Solving the Generalized Poisson Equation Using the Finite-Difference Method (FDM)". Technical Report (2012).
#   https://my.ece.utah.edu/~ece6340/LECTURES/Feb1/Nagel%202012%20-%20Solving%20the%20Generalized%20Poisson%20Equation%20using%20FDM.pdf

function V = poissonFDM(V, BC, EPS, RHO, h)

    EPS_0 = 8.854E-12 # permitivitty of free space; F/m

    # extract simulation domain size
    [Ny, Nx] = sizeof(V);

end
