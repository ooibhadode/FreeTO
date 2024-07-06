# FreeTO

 **FreeTO** stands for Freeform Topology Optimization (3D topology optimization through a structured mesh with smooth boundaries)
 and is an open-source Matlab code for initializing, optimizing, and post-processing free-form 3D topology optimization problems.

<img src="STLs/hand_TO.gif" width = "500" alt="">

 ## Syntax
 * **FreeTO(*file*, *force1*, *MeshControl*, *volfrac*, 'PropertyName', VALUE,...)** Performs topology optimization on the STL domain provided in *file* with the load-carrying STL subdomain provided in *force1*, *MeshControl* to create a structured hexahedral mesh, *volfrac* for the required volume fraction, property names and values (support conditions and force magnitudes are essential), and other default parameters. The property name-value inputs are given in Table 1.

### Table 1: Name-value inputs

--------------------------------------------------------------------------------------------------------------------
 | Name-value                        | Description                                                                                                                                  | Default value |
 |-----------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------|---------------|
 | fixed                             | Mesh file (.stl) of the regions in domain with a Dirichlet displacement boundary condition fixed in all three directions                     | -             |
 | xfixed, yfixed, zfixed            | Type of topology optimization approach                                                                                                       | -             |
 | force2 – force10                  | Mesh files (.stl) of regions in domain that carry up to 9 additional load cases                                                              | -             |
 | keepdom                           | Mesh file (.stl) of regions in domain that should remain unchanged throughout the optimization                                               | -             |
 | keep_BC                           | Specify if the load and support (fixed in all directions) regions should be included in the optimized structure                              | Yes           |
 | keep_BCx, keep_BCy, keep_BCz      | Specify if the support regions fixed in x, y, and z directions, respectively, should be included in the optimized structure                  | No            |
 | YoungsModulus                     | Specify the Young’s Modulus of the material                                                                                                  |
| PoissonsRatio                     | Specifies the Poisson’s ratio of the material                                                                                                | 0.3           |
 | optimization                      | Specify the optimization method (SIMP or SEMDOT)                                                                                             | SIMP          |
 | Fmagx, Fmagy, Fmagz               | Indicates the x, y, and z components of the loads                                                                                            | 0, 0, 0       |
 | modelName                         | Indicates the name of the mesh file (.stl) generated from the optimized structure                                                            | -             |
 | loadtype                          | Indicates if a load should be applied at a point (center node) or distributed across all nodes in a load case                                | distributed   |
 | Symmetry1, Symmetry2, Symmetry3   | Indicates a symmetry plane to mirror the optimized topology. Symmetry can be indicated a maximum of 3 times.                                 | -             |
 | direction1, direction2, direction3 | Indicates whether the mirror should occur to the left or right of the symmetry plane. It can be invoked the same number of times as Symmetry | Right         |

 ### Table 2: Outputs
 | Variable | Description                                             |                                    
 |----------|---------------------------------------------------------|
 | comp     | Compliance of optimized topology                        |
 | finalvol | 	Volume fraction of optimized topology                  |
 | eleden   | Element density array                                   |
 | gridden  | Grid density array                                      |
 | elenum1  | Number of active elements in the optimization           |
 | elenum2  | Total number of elements (active + passive)             |
 | 	        | Display of the optimized topology                       |
 | 	        | Display of the compliance and volume fraction histories |   
 | 	        | STL file (optional)                                     |   

## Examples:

The following examples show how **FreeTO** is used.
Input STL files of the sample design problems have been prepared and are in the STLs folder in this repository. The input file preparation is shown in a video for Example 1 [here](https://1drv.ms/v/s!ArlEZ0UGGKnchvZ0yWGqPRQcVgRg1Q?e=mK8vHb).

  *  Example 1:
    Optimize a GE bracket given the freebody diagram on the left of Figure 1

     - **FreeTO('STLs\GE_domain.stl','STLs\GE_force.stl',80,0.3,'fixed','STLs\GE_fixed.stl','force2','STLs\GE_force.stl','Fmagy',[0 -2000],'Fmagz',[1500 0],'YoungsModulus',210e9)**

<img src="STLs/Fig1.png" width = "500" alt="">

  *  Example 2:
     Optimize one half of an airplane bracket 

     - **FreeTO('STLs\air_domain.stl','STLs\air_force.stl',90,0.2,'fixed','STLs\air_fixed.stl','zfixed','STLs\air_zfixed.stl','force2','STLs\air_force.stl','force3','STLs\air_force.stl','Fmagx',[1000 1324 0],'Fmagy',[0 -1324 -2500],'YoungsModulus',210e9,'Symmetry1',"x-y",'direction1',"right",'optimization','SEMDOT','modelName','STLs\air_brack_TO.stl')**

  <img src="STLs/Fig2.png" width = "500" alt="">

  *  Example 3:
     Optimize a human hand model with loads on the fingertips

     - **FreeTO('STLs\hand_domain.stl','STLs\hand_force1.stl',90,0.3,'fixed','STLs\hand_fixed.stl','force2','STLs\hand_force2.stl','force3','STLs\hand_force3.stl','force4','STLs\hand_force4.stl','force5','STLs\hand_force5.stl','Fmagz',2e3*ones(1,5),'YoungsModulus',210e9)**

<img src="STLs/Fig3.png" width = "500" alt="">

  *  Example 4:
     Optimize the quarter of a quadcopter

     - **FreeTO('STLs\quad_domain.stl','STLs\quad_force1.stl',70,0.3,'fixed','STLs\quad_fixed.stl','xfixed','STLs\quad_xfixed.stl','yfixed','STLs\quad_yfixed.stl','force2','STLs\quad_force2.stl','force3','STLs\quad_force3.stl','Fmagz',[-1500 -1500 -1000],'YoungsModulus',2e9,'Symmetry1','y-z','Symmetry2','z-x','direction1','right','direction2','left','modelName','STLs\quad_TO.stl','keep_BCz','yes')**

<img src="STLs/Fig4.jpg" width = "500" alt="">


## Supporting Open-Source Codes
**FreeTO** utilizes other open-source codes such as [intriangulation.m](https://www.mathworks.com/matlabcentral/fileexchange/43381-intriangulation-vertices-faces-testp-heavytest) by Johannes Korsawe, and modifications to [top3D.m](https://www.top3d.app/) by Liu and Tovar, [SEMDOT.m](https://blogs.deakin.edu.au/dot/software/) by Fu et al., and [phi2stl.m](https://asmedigitalcollection.asme.org/computingengineering/article/17/4/041012/474348/An-Open-Source-Framework-for-Integrated-Additive) by Vogiatzis et al. FreeTO on this repository uses the optimality criteria method (OCM) but works much better with the Method of Moving Asymptotes (MMA) as optimizer. The implementation of MMA is given in SIMP.m and SEMDOT.m but the lines are commented out as the codes are not on this repository and the user should request them from [Prof. Krista Svanberg](https://people.kth.se/~krille/Welcome.html) by email. If using MMA, comment out lines 79 to 84, 135, 137 and uncomment lines 44 to 56, 86 to 99 in SEMDOT.m and SIMP.m. Furthermore, in SEMDOT.m and SIMP.m, tolx = tol_thresh = 1e-3, beta = ER = 0.5 are recommended values for MMA.

## To Cite
If you find this code helpful in your work, please cite [this open access paper](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4876754)
