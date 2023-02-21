# BO_STLS
## Workflow of the proposed identification method
### Step 1, motion data collection and pre-processing. 
In the whole modeling process, only the data of UMV speed components and rudder angle that can be easily collected by on-board sensors are required to train the model. Data on UMV movements are collected using sensors such as IMU/GPS/USBL, while their control signals are recorded through the control system. The data obtained from measurements in real environments can be contaminated. In sea trials they are disturbed by environmental factors such as wind and wave currents; even indoor pool tests are affected by ship vibrations, mechanical noise etc. The 'wavelet' toolbox of MATLAB performs well in ship maneuverability prediction .
### Step 2, divide the training and validation sets. 
To enhance the generalization of the model, 30% of the data were selected as the validation set. Hyperparameters were obtained through a Bayesian optimizer. The dictionary library for the surface vehicles is selected; while the dictionary library for the underwater vehicles is selected.
### Step 3, train and validate the model. 
The STLS algorithm is used to discover the truly active terms from the hydrodynamic dictionary library; and the identification results are compared with the parametric and non-parametric model identification results respectively to verify the effectiveness of the method.

The above process is summarized in the flow chart
