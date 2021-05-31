## Detailed Suspension Control Co-Design Problem

<p align="center">
  <img src="detailed-suspension-acc/assets/acc2021-image-2.svg">
</p>

### Key Files

##### Simultaneous Strategy

- For a single problem instance with the simultaneous CCD strategy, run ![DSuspensionSimultaneous.m](detailed-suspension-acc/simultaneous/DSuspensionSimultaneous.m)
- For many problems with different options (tolerances, number of time points, derivative methods, etc.) with the simultaneous CCD strategy, run ![DSuspensionSimultaneous_Sens.m](detailed-suspension-acc/simultaneous/DSuspensionSimultaneous_Sens.m)

##### Nested Strategy

- For a single problem instance with the nested CCD strategy, run ![DSuspensionNested.m](nested/DSuspensionNested.m)
- For many problems with different options (tolerances, number of time points, derivative methods, etc.) with the nested CCD strategy, run ![DSuspensionNested_Sens.m](nested/DSuspensionNested_Sens.m)

### Reference

##### Code and results for the following reference:

A. K. Sundarrajan and D. R. Herber, "*Towards a Fair Comparison Between the Nested and Simultaneous Control Co-Design Methods Using an Active Suspension Case Study*", American Control Conference, May 2021

##### Original reference:

J. T. Allison, T. Guo, and Z. Han, "*Co-Design of an Active Suspension Using Simultaneous Dynamic Optimization*", Journal of Mechanical Design, vol. 136, no. 8, Jun. 2014, doi: 10.1115/1.4027335. [Online]. Available: http://dx.doi.org/10.1115/1.4027335 

### Formulation
<p align="center">
  <img height="225" src="detailed-suspension-acc/assets/acc2021-image-1.svg">
</p>

![formulation](detailed-suspension-acc/assets/formulation.svg)
Note that this formulation is for the simultaneous problem. The nested strategy (inner/outer loops with fixed plant design for the inner-loop problem) is also implemented as indicated above.

<!-- ### Solution -->