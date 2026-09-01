# OpenWAM-Refactored
OpenWAM refactored is a 1D gas dynamic code primarily targeting the simulation of internal combustion engines, forked from the original [repository](https://github.com/CMT-UPV/OpenWAM) made by CMT-UPV.

## Brief Overview
The ui allows the user to place nodes representing valves, pistons, pipes, and other components, assembling them into a model that resembles to interior mapping of an internal combustion engine. Each component has individual parameters such as diameter, length, bore, stroke, ambient air temperature/pressure, vehicle weight, etc., to allow the simulation of very specific situations. Quantities like pressure, temperature, power, and others can be measured at specific positons anywhere in the graph, allowing precise results gathering for comparison to research or the real world.

## Notable Changes
- Added a 2 zone scavenging model to the two stroke path
- Replaced the traditional Roe interior solver and Benson Junctions with a RoeM interior + GJM junctions to improve accuracy and quell issues with energy conservation
- Added multithreading to the common simulation path, allowing multiple pipe interiors to be processed in parallel
- Tested and documented the majoriy of ui features
- Fixed many bugs and improved error messages

## Get Started
- Users: Download the latest release for you architecture and read the user guide.
- Developers: Clone the repository and read the user and developer guides in the tutorial directory.

## AI Usage
Claude code was employed to greatly increase the speed of isololating bugs, directly writing solutions, and generating tutorial pdfs. Research, planning, the creation of examples, and everything else was done by a human (myself).

## Citations

Original OpenWAM [website](https://openwam.webs.upv.es/docs/)

LÓPEZ, J. M. P. (2017). STUDY AND UNDERSTANDING OF THE SOFTWARE OPENWAM AND USE FOR ENGINES SIMULATION, Universitat Politècnica de València. https://www.scribd.com/document/420809839/openwam

Qiao, Y., Duan, X., Huang, K., Song, Y., & Qian, J. (2018). Scavenging Ports’ Optimal Design of a Two-Stroke Small Aeroengine Based on the Benson/Bradham Model. Energies, 11(10), 2739. https://doi.org/10.3390/en11102739

Beccari, S., Pipitone, E., & Genchi, G. (2015). Calibration of a Knock Prediction Model for the Combustion of Gasoline-LPG Mixtures in Spark Ignition Engines. Combustion Science and Technology, 187(5), 721–738. https://doi.org/10.1080/00102202.2014.960925

Choi, S., Kolodziej, C., Hoth, A., and Wallner, T., "Development and Validation of a Three Pressure Analysis (TPA) GT-Power Model of the CFR F1/F2 Engine for Estimating Cylinder Conditions," WCX World Congress Experience, Detroit, Michigan, United States, April 10, 2018, https://doi.org/10.4271/2018-01-0848.

Jan, S., Mohammed, A., Elkhazraji, A., Masurier, J., et al., "The Road Towards High Efficiency Argon SI Combustion in a CFR Engine: Cooling the Intake to Sub-Zero Temperatures," WCX SAE World Congress Experience, Detroit, Michigan, United States, April 21, 2020, https://doi.org/10.4271/2020-01-0550.

Filipi, Z., Wang, Y., and Assanis, D., "Effect of Variable Geometry Turbine (VGT) on Diesel Engine and Vehicle System Transient Response," SAE 2001 World Congress, Detroit, Michigan, United States, March 5, 2001, https://doi.org/10.4271/2001-01-1247.

Filipi, Z., Wang, Y., & Assanis, D. (2004). Variable geometry turbine (VGT) strategies for improving diesel engine in-vehicle response: a simulation study. International Journal of Heavy Vehicle Systems, 11(3/4), 303. https://doi.org/10.1504/ijhvs.2004.005453

Capata, R., & Sciubba, E. (2021). Study, Development and Prototyping of a Novel Mild Hybrid Power Train for a City Car: Design of the Turbocharger. Applied Sciences, 11(1), 234. https://doi.org/10.3390/app11010234

Capata, R. (2021). Experimental Fitting of Redesign Electrified Turbocompressor of a Novel Mild Hybrid Power Train for a City Car. Energies, 14(20), 6516. https://doi.org/10.3390/en14206516

Lombardi, S., Ricci, F., Martinelli, R., Tribioli, L., Grimaldi, C. N., & Bella, G. (2023). Energy Analysis of a Novel Turbo-Compound System for Mild Hybridization of a Gasoline Engine. Energies, 16(18), 6444. https://doi.org/10.3390/en16186444

Corberán, J. M., & Gascón, M. Ll. (1995). TVD schemes for the calculation of flow in pipes of variable cross-section. Mathematical and Computer Modelling, 21(3), 85–92. https://doi.org/10.1016/0895-7177(94)00216-b

Kim, S., Kim, C., Rho, O.-H., & Kyu Hong, S. (2003). Cures for the shock instability: Development of a shock-stable Roe scheme. Journal of Computational Physics, 185(2), 342–374. https://doi.org/10.1016/s0021-9991(02)00037-2

Hong, S. W., & Kim, C. (2011). A new finite volume method on junction coupling and boundary treatment for flow network system analyses. International Journal for Numerical Methods in Fluids, 65(6), 707–742. https://doi.org/10.1002/fld.2212

Serrano, J. R., Arnau, F. J., Piqueras, P., & García-Afonso, O. (2013). Application of the two-step Lax and Wendroff FCT and the CE-SE method to flow transport in wall-flow monoliths. International Journal of Computer Mathematics, 91(1), 71–84. https://doi.org/10.1080/00207160.2013.783206
