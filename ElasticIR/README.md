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
Say you have $`x`$, a vector of $`n`$ state variables, $`y`$, a vector of $`k`$ controls, and $`f`$, a vector of $`m=n+k`$ equations. Let $`\textbf{x}=(x,x')`$: and $`\textbf{y}=(y,y')`$

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
& 
f(\textbf{x},\textbf{y})&=\begin{bmatrix}
f_1(\textbf{x},\textbf{y}) \\
\vdots \\
f_m(\textbf{x},\textbf{y})
\end{bmatrix}
\end{aligned}
$$

Schmitt-Grohé and Uribe made a contribution to the literature by integrating Matlab's symbolic toolbox with this standard economic paradigm, allowing for the formulation of Taylor-type approximations to policy functions and state evolution equations. Specifically, we denote $` h(x) `$ as the state evolution equation and $`g(x)`$ as the policy functions. To do a second order approximation of the evolution of the $` i`$th  state variable, we do the following 
$$
h(x_i) \approx x_i^*+h_x^*(x_i)\widehat{x}+\frac{1}{2}\left[\widehat{x}^Th_{xx}^*(x_i)\widehat{x}+\sigma^2h^*_{\sigma\sigma}(x_i)\right]
$$
Notation: we have the usual $`x_i^* `$ denoting the steady state of $`x_i`$, $`\widehat{x}_i`$ denoting "linearization" by $`x_i-x_i^* `$, and $` \sigma^2 `$ denoting the variable of the stochastic component of the system. The non-standard notation is are the coefficient components. These correspond to the coefficients necessary to approximate the solution to the system of equations at the steady state. Specifically, let $`\textbf{x}^h=(x,h(x))`$: and $`\textbf{y}^g =(g(x),g(h(x)))`$. We can exploit the fact that the $`h`$ and $`g`$ must satisfy the following 

$$
\begin{aligned}
x'&=h(x)+\eta \sigma \varepsilon'
& 
y&=g(x)
& 
f(\textbf{x}^h,\textbf{y}^g)&=\begin{bmatrix}
f_1(\textbf{x}^h,\textbf{y}^g) \\
\vdots \\
f_m(\textbf{x}^h,\textbf{y}^g)
\end{bmatrix}
\end{aligned}
$$

to compute their derivatives. Two tricks in particular help simplify things: if a function $`F`$ evaluated at all values possible values of its argument is 0, then any derivative of $`F`$ also must be 0. This is the exact situation we have with our system. When we put our vector of functions $`f`$ into Matlab's symbolic toolbox, we move everything over to the RHS so the LHS is equal to 0. The other trick is the toolbox treats  $`\sigma`$ as an implicit argument to $`h`$ and $`g`$. Using this approach, we can find the derivatives of these functions and then evaluate them at the steady state to construct a Taylor approximation. These values are our coefficient components from before. 

Now we have the context to decompose these coefficient matrices. The notation itself tries to enunciate the mechanical components of this approximation. 
* To lay the syntactic background and simplify discussion on second derivatives, let's call $`h_x`$ the matrix of first derivatives of $`h`$. So each element in this matrix 
* $`h_x^*(x_i)`$ is the component relating to $`x_i`$ in the matrix of first derivatives of $`h`$ evaluated at the steady state (thus the star). To be clear: the component relating to $`x_i`$ is the $`i`$th row of the first derivative matrix, representing the partial derivative w.r.t $`x_i`$. Therefore,  $`h_x^*(x_i)`$ is a $`n `$-dimensional row vector. 
* For the second derivatives, we technically don't have a matrix but instead an array of second derivatives: a collection of $`n`$ matrices of size $`n \times n`$. 
