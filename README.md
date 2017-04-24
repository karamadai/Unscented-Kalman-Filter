----------

Unscented Kalman Filter
----------

The goal of this project is to develop an Unscented Kalman Filter (UKF) to estimate the state of an object (bicycle) detected using a Radar and Lidar sensor. The UKF filter uses a non linear process model called the "Constant Turn Rate and Velocity" (CTRV) to describe the motion of the object. 

#### Code Structure

The code reads the sensor data from the text file obj_pose-laser-radar-synthetic-input.txt  
The input file provides the Radar and the Lidar data in the following format:

L(for laser) meas_px meas_py timestamp gt_px gt_py gt_vx gt_vy
R(for radar) meas_rho meas_phi meas_rho_dot timestamp gt_px gt_py gt_vx gt_vy

The code outputs to the console the accuracy of the predictions as the root mean squared error (RMSE)  px, py, vx and vy . 

The predictions (estimations) and the NIS values corresponding to the sensor inputs are saved to the output file output.txt.

The program can be run in 3 modes by setting the values of use_radar_ and user_laser_ variables (in ufk.cpp (line 16, 18) to true or false.

#### Initialization of Covarience matrix P_ and State vector x_
The Covariance matrix P_ and the state vector is initialized as follows:

P_:  All the diagonal elements of 5x5 matrix is set to 0.85 and the rest of the elements is set to 0 (ukf.cpp, line 28)
x_: The state vector is initialized to px and py values of the first sensor reading (ufk.cpp, line 85). The rest of the elements are set to 0.

#### Process Noise tuning 
The process noise  standard deviation of the linear acceleration for the bicycle is set to 0.5 rad/sec^2
The process noise standard deviation of the yaw acceleration for the bicycle  is set to 1.2 rad/s^2
The radar and laser noise measurements are defined in lines 40 to 51 in ukp.cpp.

#### Normalized Innovation Squared (NIS) consistency check
The consistency of the UKF filter is measured and tuned using Normalized Innovation Squared (NIS) method. The UKF filter estimation follows a Chi-Squared distribution. The following figures shows the UKF estimation for Laser and Radar sensor with 2 and 3 degress of freedom respectively:

![enter image description here](https://github.com/karamadai/Unscented-Kalman-Filter/blob/master/NIS_Laser.PNG?raw=true)

![enter image description here](https://github.com/karamadai/Unscented-Kalman-Filter/blob/master/NIS_Radar.PNG?raw=true)




#### Filter Performance 
The performance of the filter is measured by the root mean squared error of the estimated state against the ground truth.  

The filter gives the following RMSE for the input data in obj_pose-laser-radar-synthetic-input.txt . The output below shows the UKF filter estimations are much more closer to the ground truth when the Radar and Lidar data is combined.
 ________________________________________________________________
| RSME Radar Only  | RSME Lidar    | RSME Combined   |
|------------------|---------------|-----------------|
| 0.15             | 0.10          |  **0.058**      |
| 0.21             | 0.09          |  **0.087**      |
| 0.27             | 0.57          |  **0.22**       |
| 0.28             | 0.25          |  **0.21**       |

#### UKF and EKF Comparision
The key difference between the UKF and EKF filter is in the way the two handles non linear processes. The EKF filter approximates the non linear function using the Taylor's polynomial series truncated at the first order partial derivative (The Jacobian) to estimate the state vector at time k+1. 
The UKF on the other hand identifies sigma points that adequately represents the co-variance distribution of the state vector at time k and estimates the values of these sigma points at time k+1 by evaluating it using the non-linear function.  The new co-variance and the new mean at k+1 are then derived from the estimated sigma at k+1.

The following table compares the RMSE values of UKF and EFK filters. The UFK shows higher accuracy in the velocity Vx and Vy

| RMSE EKF   | RMSE UKF  |
|-----------|-----------|
|  0.06     |   0.058   |
|  0.06     |   0.087   |
|  0.53     |   **0.22**    |
|  0.54     |   **0.21**    |

> Written with [StackEdit](https://stackedit.io/).
