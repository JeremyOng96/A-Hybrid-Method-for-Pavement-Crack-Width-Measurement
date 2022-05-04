# A-Hybrid-Method-for-Pavement-Crack-Width-Measurement
A method to calculate the width of a binary image as described in: https://www.sciencedirect.com/science/article/abs/pii/S026322412200505X#:~:text=The%20hybrid%20method%20obtains%20the,close%20to%20the%20orthogonal%20direction.

The method identifies a group of points that are close to the orthogonal projection vector and selects a pair of points that is the shortest as shown below:<br>
<p align='center'>
  ![comparison of measurements](/comparison_measurement.png)
 <p>
   
The benefits of the hybrid method are as follows:
   (1) Allows the user to tune the strength of the orthogonal projection and the shortest method
   (2) Prevents overestimation effects by the orthogonal projection method when cracks are curved
   (3) More distinct crack width measurements than the shortest method
   (4) More robust to binary images that have unparallel boundaries
   
How to use:
   (1) Input binary image
   (2) Tune the pruning strength
   (3) Tune the balancing coefficient

   
