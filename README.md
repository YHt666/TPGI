# TPGI
TPGI is an efficient solver for the inverse kinematics of cable-driven manipulators with pure rolling joints using a geometric iterative approach.

##  About the Author
Haotian Yang received the B.S. degree in Mechanical Design, Manufacturing and Automation, from Harbin Institute of Technology, Shenzhen (HITSZ), Shenzhen, China, in 2022. He is currently working toward the M.S. degree in Electronic Information (Artificial Intelligence) with Center of Intelligent Control and Telescience, Tsinghua Shenzhen International Graduate School, Tsinghua University, Shenzhen, China. His research interests include kinematics and dynamics of robots, robotic manipulation, and deep reinforcement learning.

## Dependencies
TPGI is implemented in MATLAB R2020a. The Robotics Toolbox for MATLAB (10.4 version) is used, which can be downloaded from https://petercorke.com/toolboxes/robotics-toolbox/.

## Usage

### TPGI solver
The source code of TPGI and a test program are provided in the *TPGI_solver* folder.

### Statistical experiments
The source code of the KDL method, the SQP method and the TPGI method is given in the *Comparison_of_multiple_IK_methods* folder. A program for performance comparison between the three algorithms is also provided. In the *Runtime_distribution_and_workspace_coverage* folder, the program for plotting the runtime distribution is given in two scripts according to the running order.

### Generality tests
The *Generality_tests* folder contains parameters for the four manipulators to be tested, corresponding TPGI solvers, and performance test programs.

## Notes
The solvers' parameters (e.g., the max runtime) need to be adjusted according to the hardware parameters of the running platform！
