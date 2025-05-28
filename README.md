# Double Driven Cavity Flow with Oscillating Lids (MATLAB)

This repository contains MATLAB scripts for simulating and analyzing the 2D flow dynamics within a square cavity where both the upper and lower lids oscillate at potentially different frequencies. The primary goal of this numerical study is to investigate how these oscillations influence fluid mixing within the cavity. This project was developed as part of the Numerical Fluid Dynamics course under Prof. G. Coppola at the University of Naples “Federico II.”

## Project Overview

The simulation solves the 2D incompressible Navier-Stokes equations in primitive variables (velocity components u, v, and pressure p) using a finite difference method on a staggered Harlow-Welch grid.

Key features and methodologies include:

*   **Governing Equations:** 2D Navier-Stokes equations for incompressible flow.
*   **Numerical Scheme:**
    *   **Spatial Discretization:** Staggered grid (Harlow-Welch) for velocity and pressure.
    *   **Projection Method:** A classic projection method (Fractional Step Method) is used to decouple pressure and velocity calculations.
*   **Time Integration:** Explicit Runge-Kutta 4th order (RK4) scheme.
*   **Boundary Conditions:**
    *   No-slip conditions on stationary vertical walls.
    *   Time-varying horizontal velocity on the upper and lower lids:
        *   Upper lid: `U_up(t) = U_ref * sin(ω_ref * t)`
        *   Lower lid: `U_dwn(t) = U_ref * cos(K * ω_ref * t)`
    *   `U_ref` is a reference velocity.
    *   `ω_ref` is the reference angular frequency (typically for the upper lid).
    *   `K` is the ratio of the lower lid's angular frequency to the upper lid's angular frequency.
*   **Particle Tracking:**
    *   The paths of nine distinct particles, initially clustered near the center of the cavity, are tracked over time using the same RK4 scheme.
    *   Velocity at particle locations is obtained via 'spline' interpolation from the staggered grid.
*   **Mixing Quantification:**
    *   An 'index of dispersion' (`I/I_0`) is calculated. This index is based on the moment of inertia of the particle system with respect to its centroid, normalized by its initial value. It provides a quantitative measure of how spread out the particles become.
*   **Parametric Studies:** The code is designed to investigate the influence of:
    *   Reynolds number (`Re`)
    *   Lid frequency ratio (`K`)
    *   Reference angular frequency (`ω_ref`)
*   **Flow Visualization:**
    *   Particle scattering plots.
    *   Streamline plots (using `streamslice`).
    *   Stream function (`Ψ`) contour maps (calculated from the velocity field).
    *   Dispersion index (`I/I_0`) vs. time plots.

## File Structure

*   `DrivenHW_oscillating_lid_ANTONIO_CAROTENUTO.m`: The main script that sets up the simulation parameters, initializes variables, runs the time-stepping loop, performs particle tracking, and generates plots.
*   `RHS_HW_OL.m`: Function to calculate the Right-Hand Side (RHS) of the momentum equations. This includes the convective and diffusive terms. It's called at each stage of the RK4 integrator.
*   `DivCalc.m`: Function to calculate the divergence of the velocity field on the staggered grid.
*   `GradCalc.m`: Function to calculate the gradient of the pressure field on the staggered grid.
*   `GivePsi.m`: Function to calculate the stream function (`Ψ`) from the velocity components (U, V) by solving a Poisson equation for vorticity.
*   `plot_disp.m`: A separate script used for post-processing and plotting the dispersion index results from saved `.mat` files for various parameter combinations.


## How to Run

1.  **Main Simulation:**
    *   Open `DrivenHW_oscillating_lid_ANTONIO_CAROTENUTO.m` in MATLAB.
    *   Adjust simulation parameters within the script as needed:
        *   `Lx`, `Ly`: Domain dimensions (currently 1x1).
        *   `Nx`, `Ny`: Number of grid points (currently 80x80).
        *   `Re`: Reynolds number.
        *   `T`: Total simulation time.
        *   `Uref`: Reference lid velocity.
        *   `Omega_ref`: Reference angular frequency (`ω_ref`).
        *   `K`: Lid frequency ratio.
    *   Run the script.
    *   Figures will be generated showing particle scattering, streamlines, stream function contours, and the dispersion index evolution.

## Simulation Details (from `DrivenHW_oscillating_lid_ANTONIO_CAROTENUTO.m`)

*   **Domain:** Square, `Lx = 1`, `Ly = 1`.
*   **Grid:** `Nx = 80`, `Ny = 80` (uniform).
*   **Initial Particle Positions:** A 3x3 cluster of 9 particles around the cavity center:
    `(Lx/2, Ly/2)`, `(Lx/2, 9*Ly/20)`, `(Lx/2, 11*Ly/20)`, etc.
*   **Default Parameters (can be changed in the script):**
    *   `Re = 2000`
    *   `T = 30` (seconds)
    *   `Uref = 1`
    *   `Omega_ref = 0.8`
    *   `K = 1`

## Notes

*   The code uses global variables for convenience in accessing parameters within functions.
*   The Poisson equation for pressure is solved using MATLAB's backslash operator (`\`) with a discrete Laplacian operator constructed using `numgrid` and `delsq`.
*   Particle positions are constrained to remain within the domain boundaries.

## Author

*   A. Carotenuto

## Acknowledgements

*   Prof. G. Coppola, for supervision and guidance during the Numerical Fluid Dynamics course.

## ToDo / Potential Improvements

*   Implement more advanced boundary condition handling (e.g., for outflow).
*   Explore adaptive time-stepping.
*   Parallelize parts of the computation for larger grids.
*   Add more comprehensive error checking and input validation.
*   Improve code commenting and documentation within scripts.
*   Containerize the environment (e.g., using Docker) for easier reproducibility if specific MATLAB toolboxes or versions become critical.
