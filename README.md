# attitudeTimeOptimalManeuvers
Time optimal rest-to-rest attitude maneuvers are solved for using the Minimum principle

Code for solving the 1-dimensional and the 3-dimensional attitude maneuvers is provided. The main scripts are in the root folder, the loop script attempts a continuation process.
The H-infinity controller can be obtained by commenting out the denominator in the omegaControl.m function as shown in the project report.

Results for a given maneuver are included in results. The format of this data is directly compatible with the display code provided at the end of the main 3-dimensional script.