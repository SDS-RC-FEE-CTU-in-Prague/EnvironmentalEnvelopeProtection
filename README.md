# EnvironmentalEnvelopeProtection
Implementation of the environmental envelope protection algorithm in Simulink with testing environment in IPG Carmaker.

This is a source code of our implementation of the environmental envelope protection algorithm, described in a paper, which is under review process. Environmental envelope protection is implemented in the core as a model predictive controller. As a framework for the optimal control problem formulation, we used [CasADi](https://github.com/casadi/casadi); as a solver for nonlinear problems, we used [IPOPT](https://github.com/coin-or/Ipopt), which utilized ma27 linear solver from [HSL package](https://licences.stfc.ac.uk/product/coin-hsl-archive). As an Implementation environment, we used Matlab/Simulink R2021a. As a test environment, we used IPG CarMaker 10.1.

This project also contains the setting of the used car and test environments in CarMaker. Video results of the experiments are published on our [YouTube channel](https://www.youtube.com/@smartdrivingsolutions3633).

The folder 'EEP in C' also contains the implementation of the same code in C/C++. To run this code, make sure that you have installed all the required libraries:
  	1) Install the hsl https://github.com/coin-or-tools/ThirdParty-HSL
	  2) Install the asl https://github.com/coin-or-tools/ThirdParty-ASL
	  3) Install the ipopt https://github.com/coin-or/Ipopt
	  4) Install the casadi https://github.com/casadi/casadi/wiki/InstallationLinux (cmake -DWITH_IPOPT=ON -DWITH_HSL=ON)
	  5) Do not forget to locate the lib folder and add it to bashrc (in my case, for the default installation, it goes to /usr/local/lib)

If you want to contact authors, please use the following emails:
_denis.efremov@fel.cvut.cz_,
_tomas.hanis@fel.cvut.cz_.

# License and Acknowledgment 
This work is licensed under the [Attribution-NonCommercial-ShareAlike 4.0 International license](https://creativecommons.org/licenses/by-nc-sa/4.0/).

This work was supported in part by Toyota Motor Europe; in part by the Grant Agency of the Czech Technical University in Prague under Grant SGS22/166/OHK3/3T/13; in part by the Slovak Research and Development Agency under Grant APVV-20-0261 and Grant APVV-21-0019; in part by the Scientific Grant Agency of the Slovak Republic under Grant VEGA 1/0490/23; and in part by the Ministry of Education, Youth and Sports of the Czech Republic, under Project 8X20037.
