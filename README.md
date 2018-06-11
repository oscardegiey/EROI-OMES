
# Quantification of the electricity storage needs in Western Europe for the energy transition: an EROI based analysis



## Getting Started

clone the git: 

In the terminal enter

```
git clone https://github.com/sramelot/memoire2018
```
	
in the appropriate folder.

The git project is now in the desired folder.

### Prerequisites

You need to install the following packages:

* Anaconda 
* Pulp
* PrettyTable

Moreover, the actual codes requires the solver CPLEX from the CPLEX Optimization Suite of IBM. 
A student license can be obtained here: https://www.ibm.com/products/ilog-cplex-optimization-studio/pricing

However, CPLEX is not mandatory for the program to run. The line 'lp.solve(CPLEX())' in the file 'lp_problem_Pulp.py' can be changed to 
'lp.solve()' to use the default solver. However this solver will show worse perfomances than CPLEX for low time steps.


## Authors

* **Thibault Martinelle**
* **Sylvain Ramelot **


## Acknowledgments

* This code was highly inspired from the work of **Gauthier Limpens** in the article "*Electricity storage needs for the energy transition : an EROI based analysis

illustrated for the case of Belgium*"

