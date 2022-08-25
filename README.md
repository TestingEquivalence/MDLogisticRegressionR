This project belongs to a published article:

Vladimir Ostrovski:
"Testing equivalence to binary generalized linear models with application to logistic regression""
Statistics & Probability Letters,
Volume 191,
2022,
109658,
ISSN 0167-7152,
https://doi.org/10.1016/j.spl.2022.109658.
(https://www.sciencedirect.com/science/article/pii/S0167715222001778)



 Please take a look into examples.R to get started!   

 This project provides the minimum distance logistic regression for the situation
 that the response is binary (independent variable is binary)
 and all covariates (all dependent variables) are categorical.

 The minimum distance logistic regression uses the Euclidean distance between
 the observed counting frequencies p(n,x) and the corresponding logistic regression probabilities p(b,x),
 where n is the number of observation, x is a value of covariates and
 b is the vector of coefficients.
 The regression coefficients b are estimated using the minimization of this distance.

 In addition, this project provides the means to demonstrate that 
 the fitted logistic model is sufficiently close to the true underlying distribution.

 Let us consider the equivalence test problem:
 H0={ distance between the model und true underlying distrubition is larger than epsilon}
 H1={ distance between the model und true underlying distrubition is smaller than epsilon},
 where epsilon is the tolerance parameter.
 If H0 can be rejected for an appropriate value of epsilon, then the logistic regression model
 is sufficiently close to the true underlying distribution. 
 This program computes the minimum tolerance parameter epsilon, for which H0 can be rejected
 and hence the equivalence can be established.
 There are 4 different methods to compute the minimum tolerance parameter epsilon at the given 
 significance level:
 - approximation by asymptotic distribution, uses also asymptotic variance
 - approximation by asymptotic distribution using bootstrapped variance
 - empirical bootstrap method
 - studentized bootstrap method, also known as bootstrap-t method
 

 
 