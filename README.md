# Data-Processor-Background-Elimination-For-Raman-And-XRD:
The Data Processor for Background Elimination is a computational tool designed to enhance the accuracy and reliability of data with huge variance. In material science, powerful techniques are widely used in various scientific and industrial fields, such as XRD and Raman Spectroscopy. However, both techniques often encounter challenges posed by background signals, especially in nanomaterial, which can obscure or distort the target signals, leading to inaccurate interpretation and quantification of results.

The Data Processor for Background Elimination addresses this issue by implementing Asymmetric Least Squares Smoothing algorithms to remove unwanted background contributions from the data. With two parameters in the algorithm, the data representation is more reproducible. The tool is specifically tailored to the unique characteristics of these techniques, offering a comprehensive solution for researchers and analysts seeking meaningful insights and references from their experimental data.
## User Interface:
![GUI](https://github.com/AntonioZelongYan/Data-Processor-Background-Elimination-For-Raman-And-XRD/assets/138164005/ce875360-3296-44f7-bcd2-e0623094f213)
## Asymmetric Least Squares Smoothing Algorithm:
For people who want to have a paper reference: https://doi.org/10.1021/ac051370e

There exist two variables to be adjusted: "p," which accounts for asymmetry, and "λ," denoting smoothness. Both of these variables require calibration to suit the specific dataset in question. In our observations, a typically favorable range for "p" is between 0.0001 and 0.001, particularly when dealing with data featuring positive peaks. As for "λ," it's commonly effective to explore values as high as possible for Raman data, although exceptions might arise. 

Also, many thanks to the answer from Rustam Guliev who provided a solution with optimized memory usage in https://stackoverflow.com/questions/29156532/python-baseline-correction-library
## More Function:
### Find peaks:
By using a implended find-peak functions in Scipy, a list of peaks can be found by an adjusted parameter, prominence. This can only help to determine the peak roughly but for finding a more accurate result, other function may be needed.
### Range setting:
Some analyses only focus on the part of the whole scan of the data. By inputting the lower and upper range of the function, data selection can be undertaken.
### Moving average:
A parameter with a range from 0 to 15 for the moving average method to remove the noise in data.
# Remember: This APP only provide a more convience template for processing data.
