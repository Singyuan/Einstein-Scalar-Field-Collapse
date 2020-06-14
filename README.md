# Einstein-Scalar-Field-Collapse
This project is about numerical simulation the evolution of Einstein scalar field in spacetime.

## Preliminary
Choptuik presented a numerical study of spherically symmetric collapse of Einstein scalar field equations. Please refer to Choptuik [1] or https://singyuan.github.io/projects/NTU/Numerical_Relativity.html . Let a 4-metric <img src="https://render.githubusercontent.com/render/math?math=g=\alpha^2dt^2%2Ba^2dr^2%2Br^2d\Omega^2">. Then, given initial scalar field, we use Einstein scalar fields to solve functions <img src="https://render.githubusercontent.com/render/math?math=a"> and <img src="https://render.githubusercontent.com/render/math?math=\alpha">.

## Introduction
Given a initial scalar field, the behavior the evolution of scalar field will be shown. In this project, the initial scalar field is chosen as <img src="https://render.githubusercontent.com/render/math?math=\psi_A(r)=Ar^2e^{-(r-5)^2}">. Hence, the parameter <img src="https://render.githubusercontent.com/render/math?math=A"> can be chosen in "exp/myParameters".

>Note: Don't choose amplitude lager than 0.00015, which will collapse the spacetime. That is, the program cannot finish. In fact, there are many parameters can be change for initial scalar field, <img src="https://render.githubusercontent.com/render/math?math=\psi_A(r)=phi0r^{pow}e^{-(r-r0)^q}">.

There are two output in folder "result", including of evolution of scalar field and central value of the lapse function <img src="https://render.githubusercontent.com/render/math?math=\alpha">. Because some of Matlab cannot support creating video, if it happens, it will output the image of scalar field every coordinate time.

>Note: If the lapse function goes to zero, then the spacetime collapse, i.e. there is a black hole.

## Quick start
1. Visit the folder "exp".
2. The parameters can be tuned in "myParameters", including of amplitude of scalar field.
3. Run the file "Demo.m".
4. Visit the folder "result".
5. The evolution of scalar field is in folder "scalarfield".
6. The central value of lapse function is in folder "data".
7. The log text including some information of such simulation is in  "data".

>Note: When the program simulates, the progress and elapse time will be presented.

## Techniques
1. Iterative Crank Nicolson method: refer to Alcubierre et al. (2000) and Teukolsky (2000). 
2. Inner boundary problem: Evans (1984).
3. Outter boundary condition: Arnowitt et al. (1983) and Gundlach (2004).
4. Summation by part: Gustanfsson (2008) and Strand (1994).
5. Diffusion term: Kreiss and Oliger (1973).

## References
1. M. W. Choptuik, Universality and Scaling in Gravitational Collapse of a Massless Scalar Field, (1993).
2. M. Alcubierre, Introduction to 3+1 Numerical Relativity (2008).