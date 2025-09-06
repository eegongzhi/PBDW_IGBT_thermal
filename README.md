# PBDW_IGBT_thermal
This project uses the Parameterized-Background Data-Weak (PBDW) method to recover the thermal field distribution of an IGBT module, based on measured data.

This project is a direct application of the following two papers, where the PBDW method is deduced and applied in various scenarios (acoustic, etc.). In this project, we study the temperature distribution on the surface of an IGBT module.
- doi:10.1002/nme.4747
- doi: 10.1137/18m116544x

**Note**: The data in the folder "./data" is incomplete. Two more data files are required, but they are too large to be uploaded. Besides, this is a project without a license. To prevent any abusive usage of the codes, please contact eegongzhi@zju.edu.cn to unlock the hidden codes and the remaining data files.

## What does the PBDW method do?
The PBDW method is recognized as a _bridge_ between _simulation data_ and _measured data_. The goal of the PBDW method is to recover the field distribution using 1) the snapshots obtained from various FEM simulations, and 2) the measured field value at some (selected) locations.

To achieve that, we first conduct a series of FEM simulations and store the calculated temperature distribution of each FEM run as a snapshot. These simulations should cover as many operational conditions as possible, e.g., under different heat sources, convection coefficients, etc. Next, we build a background space using these snapshots. Then, we select some locations for placing the sensors that generate the measured temperature at the corresponding locations. Using the locations (as well as a little bit of information on the sensors), we generate the update space. Finally, we build the PBDW formula using these two spaces. 

All of these steps above are conducted _offline_. Once we obtain the PBDW formula, we put it _online_. The sensors acquire readings (i.e., the temperature at the corresponding locations), and we put these readings into the right-hand side of the PBDW formula. By solving the formula, we can finally estimate the field distribution over the entire IGBT module. Hence, the connection between simulation data and measured data is accomplished.

In this project, greedy algorithms are used to select the optimal snapshots to construct the background field, as well as to select the optimal sensor locations to construct the update space. However, technically, such algorithms are not compulsory; as long as we can construct these 2 spaces, any algorithm would be good. Therefore, I personally believe that the specific greedy algorithm is _NOT_ a part of the PBDW method.


## Code structure


## Results


