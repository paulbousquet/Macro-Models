# Elastic Interest Rate, Schmitt-Grohé and Uribe (2003)
In the "secondO" folder, there is a second order approximation to the elastic interest rate model in "[Closing small open economy models](https://www1.columbia.edu/~mu2166/closing_jie.pdf)". All that's needed is to run eir_run. The purpose of this code is to see how well it minimizes Euler equation errors, rather than impulse responses as the paper's original code, which has a goal of producing IRF based on a first order approximation. 

Additionally, the [code](https://www1.columbia.edu/~mu2166/closing.htm) originally provided for the model is more than 20 years old and does not work in the current version of Matlab. The following codes in the "orig" folder will run successfully. 
* Run edir_model.m once and then run edir_run.m
* Files in functions.zip may be needed if not already on path

Changes made include (I'm sure the other models would run after similar changes)
* [Updated the use of the subs function](https://www.mathworks.com/matlabcentral/answers/449408-error-using-sym-subs-too-many-input-arguments-error-in-mx_model-line-176-f-subs-f-cup-cu-0), which previously allowed 4 inputs. The code has been modified so that the order of arguments is flipped if the subs function returns the same thing (where relevant)
* There is a file to manually create matricies, but it produced several errors. There were no brakcets in the lines to create the matricies and the index to begin parsing was set at 8 for all cases, when it should have only been some, removing 2 of the matrix elements.
* The steady states weren't properly defined (the "prime" variables have been manually added)

# Second Order
Say you have $`x`$, a vector of $`n`$ state variables, and $`y`$, a vector of $`k`$ controls:

$$
\begin{aligned}
x&=\begin{bmatrix}
x_1 \\
\vdots \\
x_n
\end{bmatrix}
& 
y&=\begin{bmatrix}
y_1 \\
\vdots \\
y_k
\end{bmatrix}
\end{aligned}
$$
Schmitt-Grohé and Uribe made a contribution to the literature by integrating Matlab's symbolic toolbox with this standard economic paradigm, allowing for the formulation of Taylor-type approximations to policy functions and state evolution equations. Specifically, we denote $` h(x) `$ as the state evolution equation and $`g(x)`$ as the policy functions. To do a second order approximation of the evolution of the $` i`$th  state variable, we do the following 

$$
h(x_i) \approx x_i^*+h_x(x_i)\widehat{x}+\frac{1}{2}\left[\widehat{x}^Th_{xx}(x_i)\widehat{x}+\sigma^2h_{\sigma\sigma}(x_i) \right]
$$


