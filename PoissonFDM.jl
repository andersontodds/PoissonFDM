## PoissonFDM.jl
# Todd Anderson
# Jan 24 2022
#
# Solve the generalized Poisson Equation using the finite-difference method (FDM).
# Adapted from James Nagel's 2D Poisson solver here: https://drjamesnagel.com/codes/Poisson_FDM_Solver_2D.m
# Further reading:
#   1. Nagel, J. R, "Solving the Generalized Poisson Equation Using the Finite-Difference Method (FDM)". Technical Report (2012).
#   https://my.ece.utah.edu/~ece6340/LECTURES/Feb1/Nagel%202012%20-%20Solving%20the%20Generalized%20Poisson%20Equation%20using%20FDM.pdf
#   2. Wikipedia article on the discrete Poisson Equation:
#   https://en.wikipedia.org/wiki/Discrete_Poisson_equation
# A first attempt at porting this code to Julia is provided below.  Several key parts of this code, such as the matrix inversion "\", and
# dealing with the sparse matrix, will not work yet.
#
# Dependencies:
#   SparseArrays
#   Match

V = function poissonFDM(V_0, BC, EPS, RHO, h) #argument V → V_0?

    EPS_0 = 8.854E-12 # permitivitty of free space; F/m

    # extract simulation domain size
    (Ny,Nx) = size(V_0);

    # total number of unknowns to solve for
    L = Nx*Ny;

    # instantiate the dielectric coefficient matrices
    a0 = zeros(Ny, Nx);
    a1 = zeros(Ny, Nx);
    a2 = zeros(Ny, Nx);
    a3 = zeros(Ny, Nx);
    a4 = zeros(Ny, Nx);
    
    # short-hand index notation
    X1 = 2:Nx-1;
    Y1 = 2:Ny-1;
    X2 = 1:Nx-2;
    Y2 = 1:Ny-2;

    # compute the dielectric coefficients
    a0[Y1, X1] = -1*(EPS[Y1, X1] + EPS[Y2, X1] + EPS[Y1, X2] + EPS[Y2, X2]);
    a1[Y1, X1] = 0.5*(EPS[Y1,X1] + EPS[Y2, X1]);
    a2[Y1, X1] = 0.5*(EPS[Y1,X2] + EPS[Y1, X1]);
    a3[Y1, X1] = 0.5*(EPS[Y2,X2] + EPS[Y1, X2]);
    a4[Y1, X1] = 0.5*(EPS[Y2,X1] + EPS[Y2, X2]);

    # separate coefficients into real and imaginary components
    a0r = real(a0); #real
    a1r = real(a1);
    a2r = real(a2);
    a3r = real(a3);
    a4r = real(a4);
    a0i = imag(a0); #imaginary
    a1i = imag(a1);
    a2i = imag(a2);
    a3i = imag(a3);
    a4i = imag(a4);

    # check for charge densities
    # TODO: nargin syntax doesn't exist in Julia!  Convert to optional argument, or default argument.
    #   E.g., just set h to default value
    # if nargin == 3
    #     b = zeros(Ny*Nx,1);
    # else
    #     # insert default h-value if not specified
    #     if nargin == 4
    #         h = 1;
    #     end

    #     # normalized charge matrix, from integrating charge density over a square region with size h
    #     b = -1*RHO*h^2/EPS_0;
    #     b = b[:];   # vectorize

    # end
    b = -1*RHO*h^2/EPS_0;
    b = b[:];

    # Set the four grid corners to Dirichlet boundaries to avoid confusion in the algorithm
    # (what does this mean?).  These points do not really matter anyway (ha!).
    BC[1, 1]    = 1;
    BC[1,Nx]    = 1;
    BC[Ny,Nx]   = 1;
    BC[Ny,1]    = 1;

    # Define flags for Neumann boundaries. Each case needs to be handled in its own unique way, so it
    # helps to give eachone a unique marker.
    TOP_FLAG    = -1;
    BOTTOM_FLAG = -2;
    LEFT_FLAG   = -3;
    RIGHT_FLAG  = -4;

    # Find the Neumann BCs on the edges
    Neumann_TOP     = findall(BC[1, :]  != 1);
    Neumann_BOTTOM  = findall(BC[Ny, :] != 1);
    Neumann_LEFT    = findall(BC[:, 1]  != 1);
    Neumann_RIGHT   = findall(BC[:, Nx] != 1);

    # Set flags for Neumann boundaries
    BC[1,Neumann_TOP]       .= TOP_FLAG;
    BC[Ny,Neumann_BOTTOM]   .= BOTTOM_FLAG;
    BC[Neumann_LEFT,1]      .= LEFT_FLAG;
    BC[Neumann_RIGHT,Nx]    .= RIGHT_FLAG;

    # Initialize indices and values to fill the system matrix
    NZmax = 5*L;            # maximum possible number of nonzero elements
    I   = zeros(NZmax, 1);  # i-indices
    J   = zeros(NZmax, 1);  # j-indices
    Sr  = zeros(NZmax, 1);  # real values
    Si  = zeros(NZmax, 1);  # imaginary values
    idx = 1;                # nonzero entry index

    # Begin iterating over the unknowns and prepare to fill the system matrix.
    # Filling the A-matrix is MUCH faster if you know the indices and values
    # "a priori" for ALL non-zero elements. So rather than directly fill in the
    # A-matrix, this loop stores all the nonzero (i,j) indices along with their
    # values. The A-matrix is then instantiated directly from this information.
    #   Note that by convention, MATLAB scans matrices columnwise when only a 
    # single element is indexed. So rather than wasste time vectorizing the BC-
    # matrix or initial-valued V-matrix, just use a single index to load and 
    # store values and remember this convention.
    #
    # NOTE: switch/case statements do not exist in Julia Base!  Use Match.jl instead.
    # NOTE: it looks like this is essentially a sparse array generator; perhaps try
    #   Sparse Array Constructor macro here: 
    #   https://kmsquire.github.io/Match.jl/latest/#Mathematica-Inspired-Sparse-Array-Constructor-1
    display("Filling System Matrix");
    for n in 1:L
        # Check boundary condition
        @match BC[n] begin # note n is a linear index of matrix BC
            
            # Dirichlet boundary
            1 => begin
                # A[n,n] = 1
                I[idx] = n;
                J[idx] = n;
                Sr[idx] = 1;
                idx = idx + 1;
                
                # specify right-hand side as the potential at this point
                b[n] = V_0[n];
            end

            # Top Neumann boundary
            TOP_FLAG => begin
                # A[n,n] = 1
                I[idx] = n;
                J[idx] = n;
                Sr[idx] = 1;
                idx = idx + 1;

                # A[n,n+1] = -1
                I[idx] = n;
                J[idx] = n + 1;
                Sr[idx] = -1;
                idx = idx + 1;
            end

            # Bottom Neumann boundary
            BOTTOM_FLAG => begin
                # A[n,n] = 1
                I[idx] = n;
                J[idx] = n;
                Sr[idx] = 1;
                idx = idx + 1;

                # A[n,n-1] = -1
                I[idx] = n;
                J[idx] = n - 1;
                Sr[idx] = -1;
                idx = idx + 1;
            end

            # Left Neumann boundary
            LEFT_FLAG => begin
                # A[n,n] = 1
                I[idx] = n;
                J[idx] = n;
                Sr[idx] = 1;
                idx = idx + 1;

                # A[n,n+Ny] = -1
                I[idx] = n;
                J[idx] = n + Ny;
                Sr[idx] = -1;
                idx = idx + 1;
            end

            # Right Neumann boundary
            RIGHT_FLAG => begin
                # A[n,n] = 1
                I[idx] = n;
                J[idx] = n;
                Sr[idx] = 1;
                idx = idx + 1;

                # A[n,n-Ny] = -1
                I[idx] = n;
                J[idx] = n - Ny;
                Sr[idx] = -1;
                idx = idx + 1;
            end

            # Regular point; apply 5-point star
            _ => begin
                # 5-point star convention:
                # 
                #       V2          |
                #   V3  V0  V1      | Indexing direction
                #       V4          V
                #
                # REMINDER: Single-valued indexing of a matrix will scan column-wise!

                # V0 (center) term: A[n,n] = a0[n]
                I[idx] = n;
                J[idx] = n;
                Sr[idx] = a0r[n];   # real
                Si[idx] = a0i[n];   # imaginary
                idx = idx + 1;

                # V1 (right) term: A[n,n+Ny] = a1[n]
                I[idx] = n;
                J[idx] = n + Ny;
                Sr[idx] = a1r[n];   # real
                Si[idx] = a1i[n];   # imaginary
                idx = idx + 1;

                # V2 (top) term: A[n,n+1] = a2[n]
                I[idx] = n;
                J[idx] = n + 1;
                Sr[idx] = a2r[n];   # real
                Si[idx] = a2i[n];   # imaginary
                idx = idx + 1;

                # V3 (left) term: A[n,n-Ny] = a3[n]
                I[idx] = n;
                J[idx] = n - Ny;
                Sr[idx] = a3r[n];   # real
                Si[idx] = a3i[n];   # imaginary
                idx = idx + 1;

                # V4 (bottom) term: A[n,n-1] = a4[n]
                I[idx] = n;
                J[idx] = n - 1;
                Sr[idx] = a4r[n];   # real
                Si[idx] = a4i[n];   # imaginary
                idx = idx + 1;
            end

        end

    end
    # Throw out the leftover zeros
    # NOTE: this type of indexing might not work in Julia
    I = I[1:idx-1];
    J = J[1:idx-1];
    Sr = Sr[1:idx-1];
    Si = Si[1:idx-1];

    # Combine real and imaginary coefficients into one vector
    S = Sr + Si*im;

    # Fill the system matrix
    # NOTE: sparse() is a MATLAB function!  Find equivalent Julia function or macro.
    A = sparse(I,J,S,L,L,NZmax);

    # Clear out memory before performing inversion. The matrix inversion process
    # is HUGELY memory intensive, and will greedily hog up all the RAM it can get.
    # Clearing out the major variables from the workspace will at least give a few
    # extra MB of memory to this process. It won't be much, but every little bit
    # helps if we're inverting a very large system.
    #clear I J S BC a0 a1 a2 a3 a4 a0r a0i Sr Si;
    #clear a1r a1i a2r a2i a3r a3i a4r a4i Neumann_TOP Neumann_BOTTOM;
    #clear Neumann_LEFT Neumann_RIGHT EPS RHO;
    
    # Invert the matrix. The slash operator is fantastic at doing this quickly
    # for sparse matrices.
    # NOTE: take a look at SparseArrays, StaticArrays for speed in Julia version.
    display("Solving system matrix");
    V = A\b; #not a valid Julia operation!
    display("Done!");

    # Put the potentials back in the matrix
    V = reshape(V, Ny, Nx);

end

## test function
using Match
using SparseArrays

V = ones(10,10)
BC = zeros(size(V))
EPS = ones(size(V)).*EPS_0
RHO = zeros(size(V))
RHO[3,3] = 10
RHO[8,8] = -10
h = 1

V_out = poissonFDM(V,BC,EPS,RHO,h)


## test individual lines

V = ones(10,10)
BC = zeros(size(V))
EPS = ones(size(V)).*EPS_0
RHO = zeros(size(V))
RHO[3,3] = 10
RHO[8,8] = -10
h = 1

EPS_0 = 8.854E-12 # permitivitty of free space; F/m

# extract simulation domain size
(Ny,Nx) = size(V);

# total number of unknowns to solve for
L = Nx*Ny;

# instantiate the dielectric coefficient matrices
a0 = zeros(Ny, Nx);
a1 = zeros(Ny, Nx);
a2 = zeros(Ny, Nx);
a3 = zeros(Ny, Nx);
a4 = zeros(Ny, Nx);

# short-hand index notation
X1 = 2:Nx-1;
Y1 = 2:Ny-1;
X2 = 1:Nx-2;
Y2 = 1:Ny-2;

# compute the dielectric coefficients
a0[Y1, X1] = -1*(EPS[Y1, X1] + EPS[Y2, X1] + EPS[Y1, X2] + EPS[Y2, X2]);
a1[Y1, X1] = 0.5*(EPS[Y1,X1] + EPS[Y2, X1]);
a2[Y1, X1] = 0.5*(EPS[Y1,X2] + EPS[Y1, X1]);
a3[Y1, X1] = 0.5*(EPS[Y2,X2] + EPS[Y1, X2]);
a4[Y1, X1] = 0.5*(EPS[Y2,X1] + EPS[Y2, X2]);

# separate coefficients into real and imaginary components
a0r = real(a0); #real
a1r = real(a1);
a2r = real(a2);
a3r = real(a3);
a4r = real(a4);
a0i = imag(a0); #imaginary
a1i = imag(a1);
a2i = imag(a2);
a3i = imag(a3);
a4i = imag(a4);

# check for charge densities
# TODO: nargin syntax doesn't exist in Julia!  Convert to optional argument, or default argument.
#   E.g., just set h to default value
#=if nargin == 3    # RHO and h not specified → b = 0 everywhere
    b = zeros(Ny*Nx,1);
else
    # insert default h-value if not specified
    if nargin == 4  # RHO specified, but not h → set h to 1.0m
        h = 1;
    end             # other cases: nargin < 3 → can't run this function! nargin == 5 → normalize charge matrix

    # normalized charge matrix, from integrating charge density over a square region with size h
    b = -1*RHO*h^2/EPS_0;
    b = b[:];   # vectorize

end=#

b = -1*RHO*h^2/EPS_0;
b = b[:];

# Set the four grid corners to Dirichlet boundaries to avoid confusion in the algorithm
# (what does this mean?).  These points do not really matter anyway (ha!).
BC[1, 1]    = 1;
BC[1,Nx]    = 1;
BC[Ny,Nx]   = 1;
BC[Ny,1]    = 1;

# Define flags for Neumann boundaries. Each case needs to be handled in its own unique way, so it
# helps to give eachone a unique marker.
TOP_FLAG    = -1;
BOTTOM_FLAG = -2;
LEFT_FLAG   = -3;
RIGHT_FLAG  = -4;

# Find the Neumann BCs on the edges
Neumann_TOP     = findall(BC[1, :]  != 1);
Neumann_BOTTOM  = findall(BC[Ny, :] != 1);
Neumann_LEFT    = findall(BC[:, 1]  != 1);
Neumann_RIGHT   = findall(BC[:, Nx] != 1);

# Set flags for Neumann boundaries
BC[1,Neumann_TOP]       .= TOP_FLAG;
BC[Ny,Neumann_BOTTOM]   .= BOTTOM_FLAG;
BC[Neumann_LEFT,1]      .= LEFT_FLAG;
BC[Neumann_RIGHT,Nx]    .= RIGHT_FLAG;

# Initialize indices and values to fill the system matrix
NZmax = 5*L;            # maximum possible number of nonzero elements
I   = zeros(NZmax, 1);  # i-indices
J   = zeros(NZmax, 1);  # j-indices
Sr  = zeros(NZmax, 1);  # real values
Si  = zeros(NZmax, 1);  # imaginary values
idx = 1;                # nonzero entry index

# Begin iterating over the unknowns and prepare to fill the system matrix.
# Filling the A-matrix is MUCH faster if you know the indices and values
# "a priori" for ALL non-zero elements. So rather than directly fill in the
# A-matrix, this loop stores all the nonzero (i,j) indices along with their
# values. The A-matrix is then instantiated directly from this information.
#   Note that by convention, MATLAB scans matrices columnwise when only a 
# single element is indexed. So rather than wasste time vectorizing the BC-
# matrix or initial-valued V-matrix, just use a single index to load and 
# store values and remember this convention.
#
# NOTE: switch/case statements do not exist in Julia Base!  Use Match.jl instead.
# NOTE: it looks like this is essentially a sparse array generator; perhaps try
#   Sparse Array Constructor macro here: 
#   https://kmsquire.github.io/Match.jl/latest/#Mathematica-Inspired-Sparse-Array-Constructor-1
display("Filling System Matrix");
for n in 1:L
    # Check boundary condition
    @match BC[n] begin # note n is a linear index of matrix BC
        
        # Dirichlet boundary
        1 => begin
            # A[n,n] = 1
            I[idx] = n;
            J[idx] = n;
            Sr[idx] = 1;
            idx = idx + 1;
            
            # specify right-hand side as the potential at this point
            b[n] = V[n];
        end

        # Top Neumann boundary
        TOP_FLAG => begin
            # A[n,n] = 1
            I[idx] = n;
            J[idx] = n;
            Sr[idx] = 1;
            idx = idx + 1;

            # A[n,n+1] = -1
            I[idx] = n;
            J[idx] = n + 1;
            Sr[idx] = -1;
            idx = idx + 1;
        end

        # Bottom Neumann boundary
        BOTTOM_FLAG => begin
            # A[n,n] = 1
            I[idx] = n;
            J[idx] = n;
            Sr[idx] = 1;
            idx = idx + 1;

            # A[n,n-1] = -1
            I[idx] = n;
            J[idx] = n - 1;
            Sr[idx] = -1;
            idx = idx + 1;
        end

        # Left Neumann boundary
        LEFT_FLAG => begin
            # A[n,n] = 1
            I[idx] = n;
            J[idx] = n;
            Sr[idx] = 1;
            idx = idx + 1;

            # A[n,n+Ny] = -1
            I[idx] = n;
            J[idx] = n + Ny;
            Sr[idx] = -1;
            idx = idx + 1;
        end

        # Right Neumann boundary
        RIGHT_FLAG => begin
            # A[n,n] = 1
            I[idx] = n;
            J[idx] = n;
            Sr[idx] = 1;
            idx = idx + 1;

            # A[n,n-Ny] = -1
            I[idx] = n;
            J[idx] = n - Ny;
            Sr[idx] = -1;
            idx = idx + 1;
        end

        # Regular point; apply 5-point star
        _ => begin
            # 5-point star convention:
            # 
            #       V2          |
            #   V3  V0  V1      | Indexing direction
            #       V4          V
            #
            # REMINDER: Single-valued indexing of a matrix will scan column-wise!

            # V0 (center) term: A[n,n] = a0[n]
            I[idx] = n;
            J[idx] = n;
            Sr[idx] = a0r[n];   # real
            Si[idx] = a0i[n];   # imaginary
            idx = idx + 1;

            # V1 (right) term: A[n,n+Ny] = a1[n]
            I[idx] = n;
            J[idx] = n + Ny;
            Sr[idx] = a1r[n];   # real
            Si[idx] = a1i[n];   # imaginary
            idx = idx + 1;

            # V2 (top) term: A[n,n+1] = a2[n]
            I[idx] = n;
            J[idx] = n + 1;
            Sr[idx] = a2r[n];   # real
            Si[idx] = a2i[n];   # imaginary
            idx = idx + 1;

            # V3 (left) term: A[n,n-Ny] = a3[n]
            I[idx] = n;
            J[idx] = n - Ny;
            Sr[idx] = a3r[n];   # real
            Si[idx] = a3i[n];   # imaginary
            idx = idx + 1;

            # V4 (bottom) term: A[n,n-1] = a4[n]
            I[idx] = n;
            J[idx] = n - 1;
            Sr[idx] = a4r[n];   # real
            Si[idx] = a4i[n];   # imaginary
            idx = idx + 1;
        end

    end

end
# Throw out the leftover zeros
I = I[1:idx-1];
J = J[1:idx-1];
Sr = Sr[1:idx-1];
Si = Si[1:idx-1];

# Combine real and imaginary coefficients into one vector
S = Sr + Si*im;

# Fill the system matrix
# NOTE: sparse() requires SparseArrays.
A = sparse(I,J,S,L,L,NZmax);

# Clear out memory before performing inversion. The matrix inversion process
# is HUGELY memory intensive, and will greedily hog up all the RAM it can get.
# Clearing out the major variables from the workspace will at least give a few
# extra MB of memory to this process. It won't be much, but every little bit
# helps if we're inverting a very large system.
#clear I J S BC a0 a1 a2 a3 a4 a0r a0i Sr Si;
#clear a1r a1i a2r a2i a3r a3i a4r a4i Neumann_TOP Neumann_BOTTOM;
#clear Neumann_LEFT Neumann_RIGHT EPS RHO;

# Invert the matrix. The slash operator is fantastic at doing this quickly
# for sparse matrices.
# NOTE: take a look at SparseArrays, StaticArrays for speed in Julia version.
display("Solving system matrix");
V_out_test = A\b;
display("Done!");

# Put the potentials back in the matrix
V_out_test = reshape(V_out_test, Ny, Nx);
V_out_test