%% Presentation of the project

This project has been realised during an internship made under the supervision of Larisa Beilina at Chalmers University of Technology. Please cite this repository in case of usage. If you have a question or a comment, please contact me at lorentz.karmann@ens-paris-saclay.fr.

In this project, we deal with the acoustic coefficient inverse problem. We consider two situations for this problem. The first one is a reconstruction of the acoustic coefficient using internal measurements and the second one using boundary measurements.

The explanation of the situation, the methods of resolution, the results and the analyses can be found in the PDF file "Report".


%% Documentation of the repository

The algorithms are written in Matlab language. We use the library Gypsilab for the implementation of the Finite Element Methods, which has been created by Matthieu Aussal and is available at https://github.com/matthieuaussal/gypsilab. The useful files are also available in the folder gypsilab-master/gypsilab-master.

Make sure to adapt the path for each file according to the way it is downloaded on your computer.


%% gypsilab-master

Contains the library of Gypsilab made by Matthieu Aussal.
Here is the link of his Github repository:
https://github.com/matthieuaussal/gypsilab


%% Functions

Contains functions used to implement various methods described in "Report".


%% "Simulation_internal"

Contains some applications of the reconstruction using internal measurements. (see Section 2)


%% "Simulation_boundary"

Contains some applications of the reconstruction using boundary measurements. (see Section 3)


%% "Simulation_Method_1"

Contains some applications of the reconstruction using internal measurements with Method 1. (see Appendix D)


%% Simulation_1

Contains the analyses of the parameters and the methods for the reconstruction using internal measurements. (see Section 2 and Appendix E)


%% Simulation_2

Contains the analyses of the parameters and the methods for the reconstruction using internal measurements concerning the resolution of the forward problem. (see Section 2 and Appendix E)


%% Blind_inverse_solver

Contains the analyses of the reconstruction using internal measurements without the knowledge of the pseudo-frequency. (see Section 2 and Appendix F)


%% Noise

Contains the analyses of the reconstruction using noised internal measurements. (see Section 2 and Appendix G)


%% Analysis_gradient

Contains the analyses of the reconstruction using boundary measurements. (see Section 3 and Appendix J)
