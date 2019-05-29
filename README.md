# Study of Renewable Electricity Systems for Western Europe: an EROI-based analysis
This directory contains the `python3` code for optimising the production capacity mix of Western Europe by minimising the net global **EROI**. 

This code is adapted from **Sylvain Ramelot** and **Thibault Martinelle**'s Master thesis available at 

- [https://github.com/sramelot/memoire2018](https://github.com/sramelot/memoire2018) (code)
- [http://hdl.handle.net/2078.1/thesis:14328](http://hdl.handle.net/2078.1/thesis:14328) (master thesis). 

and initially inspired from the work of **Gauthier Limpens** in the article "*Electricity storage needs for the energy transition : an EROI based analysis illustrated for the case of Belgium*" available at 

- [https://doi.org/10.1016/j.energy.2018.03.180](https://doi.org/10.1016/j.energy.2018.03.180)
	
## Getting Started

Clone the git repository: 
```
git clone https://github.com/acouplet/masterthesis
```

Choose an input dataset from [this location](https://uclouvain-my.sharepoint.com/:f:/g/personal/coupleta_oasis_uclouvain_be/EhJSfKLaRD9Du2zck0EnTb4Be8fQSux3UB3F9s4qM2vnaw?e=E3klrm) and decompress the `.zip` it inside the EROI directory:  `masterthesis/EROI`.

`python3` needs to be installed on your computer. The packages needed to run the code are available in `requirements.txt`. Install  the packages using `pip3`:

```
pip3 install -r requirements.txt
```

## Configuring a test case



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
* **Sylvain Ramelot**


## Acknowledgments

* This code was highly inspired from the work of **Gauthier Limpens** in the article "*Electricity storage needs for the energy transition : an EROI based analysis
illustrated for the case of Belgium*"

