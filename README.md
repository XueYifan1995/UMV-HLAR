# UMV HLAR
This repository includes a novel, open-source, data-driven method named Hydrodynamic dictionary Library-inspired Adaptive Regression (HLAR) for unmanned marine vehicle (UMV) dynamics model identification that incorporates physical a priori knowledge; and a generic program is written in MATLAB for the method.
This method is inspired by SINDY (https://www.pnas.org/doi/10.1073/pnas.1906995116)(https://pysindy.readthedocs.io/en/latest/).
## Requirement
A computer that can run matlab software version 2020 and above.

## Feedback, bug reports, contributions
If you find this package helpful, giving a "star" to this repositry will be a happy feedback for me! If you find a bug, or have more broader kind of quession about dynamic modeling of ,please post that in the issue page. I will try hard to respond to questions.

## Usage
See 'self_library.m' to reproduce the results for unmanned surface vehicles. See 'SINDY_npsAUV_optim.m' to reproduce the results for unmanned underwater vehicles.

## Workflow of the proposed identification method
The framework of the proposed approach is as follows.
### Step 1, motion data collection and pre-processing. 
In the whole modeling process, only the data of UMV speed components and rudder angle that can be easily collected by on-board sensors are required to train the model. Data on UMV movements are collected using sensors such as IMU/GPS/USBL, while their control signals are recorded through the control system. The data obtained from measurements in real environments can be contaminated. In sea trials they are disturbed by environmental factors such as wind and wave currents; even indoor pool tests are affected by ship vibrations, mechanical noise etc. The 'wavelet' toolbox of MATLAB performs well in ship maneuverability prediction .
### Step 2, divide the training and validation sets. 
To enhance the generalization of the model, 30% of the data were selected as the validation set. Hyperparameters were obtained through a Bayesian optimizer. The dictionary library for the surface vehicles is selected; while the dictionary library for the underwater vehicles is selected.
### Step 3, train and validate the model. 
The STLS algorithm is used to discover the truly active terms from the hydrodynamic dictionary library; and the identification results are compared with the parametric and non-parametric model identification results respectively to verify the effectiveness of the method.

The above process is summarized in the flow chart
![image](https://github.com/XueYifan1995/BO_STLS/blob/master/Flow_chart.png)



