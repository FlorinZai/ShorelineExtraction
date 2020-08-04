**ShorelineExtraction<br/>**
see more at my [home page](https://weilungao.github.io)<br/>

# Abstract
The [Matlab code](https://github.com/weilungao/ShorelineExtraction) was developed for extracting the shoreline in our JGR-ES paper [(Gao et al., 2018)](https://doi.org/10.1029/2017JF004584) and GRL paper [(Gao et al., 2019)](https://doi.org/10.1029/2018GL080447), using the elevation results of Delft3D simulations. The procedures basically follow those proposed in [Shaw et al. (2008)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GL033963).<br/>

# Example
## The numerical results of D3D
I exported the numerical results of D3D in the Quickplot module of D3D. To avoid exporting the data manually, which takes a lot a time, I used a **[Matlab script](https://github.com/weilungao/ShorelineExtraction/blob/master/run_qp.m)** to export the data automatically. Then I split the "bulk" data (3 dimensions, 2 dimensions in space and 1 in time) into every single time step use a another [code](https://github.com/weilungao/ShorelineExtraction/blob/master/extract.m). Here is an example of the split [data](https://github.com/weilungao/ShorelineExtraction/blob/master/h396.mat) (2 dimensions in space), which I used for the further procedures.<br/>
<br/>Note:
* Quickplot of D3D has lot of functions which I did not know before, however, I got familiar with thoes fuctions by firstly recording the Matlab code in the GUI of Quickplot.<br/>

## Extracting the shoreline
Then I used the [Matlab code](https://github.com/weilungao/ShorelineExtraction) to extract the shoreline:<br/>

1. Put the Matlab code and the numerical results of D3D in the same path<br/>
2. Run the Matlab function *[shoreline](https://github.com/weilungao/ShorelineExtraction/blob/master/shoreline.m)* to get the shoreline<br/>

Note:
* The notes in the [Matlab code](https://github.com/weilungao/ShorelineExtraction) are currrently in Chinese, I will try to update.

## Results
The extracted results are shown below: <br/>
![Extracted shoreline](https://weilungao.github.io/assets/images/shoreline/fig23.jpg)<br/>
**Figure 1.** Procedures of the extraction of delta shoreline: (a) binary image of the bed level result where the white portion indicates bed level greater than 0 m, (b) minimum convex polygon of the cells with bed level greater than 0 m, (c) the potential shoreline grids marked as blue points, and (d) the extracted shoreline (yellow line).

# References:<br/> 
1. [Gao, W., D. Shao, Z. B. Wang, W. Nardin, P. Rajput, W. Yang, T. Sun and B. Cui (2019). 
Long-Term Cumulative Effects of Intra-Annual Variability of Unsteady River Discharge on the Progradation of Delta Lobes: 
A Modeling Perspective. Journal of Geophysical Research: Earth Surface, 124(4): 960-973.](https://doi.org/10.1029/2017JF004584) <br/>
2. [Gao, W., D. Shao, Z. B. Wang, W. Nardin, W. Yang, T. Sun and B. Cui (2018). Combined Effects of Unsteady River Discharges and Wave Conditions on River Mouth Bar Morphodynamics. Geophysical Research Letters, 45(23): 12,903-12,911.](https://doi.org/10.1029/2018GL080447) <br/>
3. [Shaw, J. B., Wolinsky, M. A., Paola, C., & Voller, V. R. (2008). An image-based method for shoreline mapping on complex coasts. 
Geophysical Research Letters, 35(12), L12405.](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008GL033963)
