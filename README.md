# 2D Numerical Simulation of Rayleigh-Bénard Convection in Enclosed Cavities

A 2D finite-difference numerical solver written in MATLAB to simulate natural convective fluid flow and heat transfer across wide ranges of Rayleigh Numbers ($1700 \le Ra \le 20000$) and Aspect Ratios ($1 \le AR \le 100$).

## 📌 Abstract & Key Findings
Convection plays a critical role in thermal systems, geophysics, and energy storage modeling. This study investigates natural convection driven by a horizontal temperature gradient in a 2D vertical enclosure with an isothermal right hot wall, isothermal left cold wall, and adiabatic top/bottom boundaries.

### Core Highlights:
* **Onset of Convection:** Flow initiation occurs near critical $Ra \approx 370$.
* **Flow Regimes:** At lower $Ra$, fluid motion is restricted near the walls with a conduction-dominated core. As $Ra$ increases, convection encompasses the entire cavity.
* **Flow Instabilities & Secondary Cells:** At $Ra = 6500$ and beyond critical aspect ratio ($AR_{cr} \cong 15$), flow instabilities set in with the formation of secondary cells.
* **Turbulence & Scaling:** High aspect ratios shrink the laminar/transition regimes, accelerating chaotic behavior at lower $Ra$.

## 📐 Numerical Methodology
* **Governing Equations:** Navier-Stokes and Energy equations for incompressible Boussinesq fluids.
* **Coupling & Discretization:**
  * **Artificial Compressibility Method:** Pressure-velocity coupling using compressibility parameter ($\beta = \sqrt{0.6}$).
  * **FTCS Scheme:** Forward Time Central Space finite-difference discretization for explicit time marching.
* **Boundary Conditions:**
  * **Velocities ($U, V$):** No-slip condition ($U = V = 0$) on all boundaries.
  * **Temperature ($T$):** Isothermal vertical walls ($T_{left} = -0.5, T_{right} = 0.5$), Neumann insulation ($\partial T/\partial y = 0$) on horizontal boundaries.
* **Convergence Criterion:** Iterates until maximum temperature residual drops below $10^{-8}$.

## 📊 Results & Visualization

| Parameter | Range |
| :--- | :--- |
| **Rayleigh Number ($Ra$)** | $1700 \le Ra \le 20000$ |
| **Aspect Ratio ($AR$)** | $1 \le AR \le 100$ |
| **Prandtl Number ($Pr$)** | $0.71$ (Air) |

### 1. Scaling & Heat Transfer Analysis
Overall heat transfer behavior across aspect ratios ($1 \le AR \le 100$):

![log Ra vs log Nu](log%20Ra%20vs%20log%20Nu.jpg)

### 2. Flow Dynamics & Secondary Cell Formation
Velocity and temperature fields illustrating secondary cells and instabilities at higher $Ra$ and $AR$:

| Velocity Fields | Velocity & Temperature Contours |
| :---: | :---: |
| ![Velocity Vector Plot](Velocity%20Vector%20Plot%20AR%2013,%20Ra%206000.jpg) <br> *Vector Field (AR = 13, Ra = 6000)* | ![Temperature Contours](Temp%20Contours%20AR%2020%20Ra%2014000.jpg) <br> *Temperature Contours (AR = 20, Ra = 14000)* |
| ![U Velocity Contours](U%20Velocity%20Contours%20AR%2020,%20Ra%206500.jpg) <br> *U Velocity Contours (AR = 20, Ra = 6500)* | ![V Velocity Contours](V%20Velovity%20Contours%20AR%2020%20Ra%208000.jpg) <br> *V Velocity Contours (AR = 20, Ra = 8000)* |

### 3. Local Nusselt Number Distribution
![Local Nu Number](Local%20Nu%20Number%20AR%2020%20Ra%208000.jpg)  
*Local Nusselt number distribution along height for AR = 20, Ra = 8000.*

## 🛠 Code & Simulation Overview

This repository is published as a portfolio showcase of numerical modeling work for MS thesis research. 

* **Language:** MATLAB
* **Primary Script:** `natural_convective_flow.m`
* **Outputs:** Computes non-dimensional velocity fields, temperature distributions, vorticity, and local/average Nusselt numbers ($Nu$).

## 📄 License & Rights

© Wajeeha Siddiqui. All rights reserved.  
This repository and its contents are for portfolio and academic viewing purposes only. No permission is granted to reproduce, distribute, or run this code without explicit consent from the author.
