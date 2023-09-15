# Software for deep nonnegative matrix factorization with the beta-divergence:
This repository contains the software for deep KL-NMF with and without minimum-volume regularization.

Copyright (c) 2023 Valentin Leplat, Nicolas Gillis, Akwum Onwunta and Le Thi Khanh Hien  <br>
Contact: ```v.leplat@skoltech.ru``` or ```Nicolas.GILLIS@umons.ac.be```


This MATLAB software reproduces the results from the following paper:

```
@unpublished{leplat-xxxxx,
  TITLE = {{Deep Nonnegative Matrix Factorization with Beta Divergences}},
  AUTHOR = {Leplat, Valentin and Le Thi Khanh, Hien and Gillis, Nicolas and Onwunta, Akwum },
  URL = {https://},
  NOTE = {preprint},
  YEAR = {2023},
  MONTH = Sept,
  KEYWORDS = {deep nonnegative matrix factorization ; $\beta$-divergences ; minimum-volume regularization ; identifiability ; multi-block nonconvex optimization ; facial feature extraction ; topic modeling ; hyperspectral unmixing},
}
```

## Acknowledgements

The baseline algorithms used in the manuscript are courtesy of their respective authors.


## Content
 
 - /Libraries : contains helpful libaries; in particular Libraries/nmfbook-master/ contains the code of the NMF toolbox from https://gitlab.com/ngillis/nmfbook/ 
 
 - /demos : contains demo files that produce tables and figures.

 - main.m : codes that allow to run desired demos.
 
 - /Utils : contains helpful files and MatLab routines to run the demos.


## Demo file
 
 A demo is available. To proceed, please type "1", "2" or "3" when running the ```main.m``` file.
 
 ## Tuning the parameters
 
 There are several parameters that you can choose:
 - The ranks r's of the deep decomposition
 - The number of iterations for the init stage and the proposed Algorithms 1 and 2
 - The accuracy for the computation of the optimal Lagrangian multipliers for MU to respect normalization constraints,
 - Parameters for ADMM procedure from Algorithm 2: rho, threshold and maximum number of iterations. Here, we allow the user to consider a new accelerated version of the ADMM procedure, see flag "options.accADMM"
 - Parameters for min-volume: \alpha_tilde helps to setup the weights \alpha_l  for the min-vol regualirzations terms. \alpha_tilde vector corresponds to the initial and desired ratio between the min-vol terms w.r.t. to the general objective function.

 
For benchmarked approaches, the parameters have been tuned according to the original works.
 
  
  ## Reproduce figures and results from the paper
  
  To do so, you need to run the ```main.m``` file. Here, a menu is available and allows you to choose which figure or table you want to generate. Each number in the table below corresponds to a set of figures.

| Number | Content                                                         |
|--------|-----------------------------------------------------------------|
| 1      | Demo for CBCL                                                   |
| 2      | Demo for topic modeling                                         |
| 3      | Benchmarking for HSI tests                                      |
