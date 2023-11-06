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
h_i(x) \approx x_i^\ast+h_x^*(x_i)\widehat{x}+\frac{1}{2}\left[ \widehat{x}^Th_{xx}^\ast(x_i)\widehat{x}+\sigma^2h^\ast_{\sigma\sigma}(x_i) \right]
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

Now we have the context to decompose these coefficient matrices. The notation itself tries to enunciate the mechanical components of this approximation. We will break this up into several pieces to lay the syntactic background and to simplify the discussion of second derivatives (which can get complicated trying to describe in words)

* To be explicit, $`h`$ is a vector of state evolution equations where $`i`$th element corresponds to $`x' = h_i(x)`$, the evolution of the $`i`$th state variable, or

  

$$
\begin{aligned}
h_x(x)&=\begin{bmatrix}
h_1(x) \\
\vdots\\
h_n(x)
\end{bmatrix}
=\begin{bmatrix}
x_1' \\
\vdots\\
x_n'
\end{bmatrix}
\end{aligned}
$$

  

* Denote $`h_x`$ the matrix of first derivatives of $`h`$. So each element in this matrix corresponds to

  

$$
\begin{aligned}
h_x&=\begin{bmatrix}
\frac{\partial h_1}{\partial x_1} & \dots & \frac{\partial h_1}{\partial x_n} \\
\vdots &\ddots & \vdots \\
\frac{\partial h_n}{\partial x_1} & \dots & \frac{\partial h_n}{\partial x_n}
\end{bmatrix}
=\begin{bmatrix}
\nabla h_1^T \\
\vdots \\
\nabla h_n^T
\end{bmatrix}
\end{aligned}
$$

  

* From our original second order equation[^1], $`h_x^*(x_i)`$ is the $`i`$th row of  $`h_x`$ evaluated at the steady state (thus the star).  As seen above, this is all possible (first) partial derivatives of $`h_i`$. 
$$
\begin{aligned}
h_x^\ast(x_i)&=[h_x^\ast]_{x_i}=\begin{bmatrix}
\frac{\partial h_i}{\partial x_1}(x^\ast) & \cdots & \frac{\partial h_i}{\partial x_n}(x^\ast)
\end{bmatrix}
\end{aligned}
$$

* For the second derivatives, we technically don't have a matrix but instead an array of second derivatives: a collection of $`n`$ matrices of size $`n \times n`$. The first matrix is the partial derivative of the $`h_x`$ matrix w.r.t $`x_1`$, the second matrix are partials w.r.t $`x_2`$ and so on, We are interested in the partial derivatives of $`\nabla h_i`$, and these will just correspond to the $`i`$th row in each matrix. Casting this in linear algebra format  (selecting these rows and then using a row combination)

$$
\begin{aligned}
h_{xx}^\ast(x_i)&=\begin{bmatrix}
\frac{\partial ^2h_i}{\partial x_1^2}(x^\ast) & \dots & \frac{\partial^2h_i}{\partial x_1 \partial x_n} (x^\ast) \\
\vdots &\ddots & \vdots \\
\frac{\partial^2 h_i}{\partial x_n \partial x_1 }(x^\ast)  & \dots & \frac{\partial^2 h_i}{\partial x_n^2}(x^\ast) 
\end{bmatrix}
=\begin{bmatrix}
 \left(\frac{\partial}{\partial x_1}\nabla h_i^T\right)(x^\ast) \\
\vdots \\
 \left(\frac{\partial}{\partial x_i}\nabla h_i^T\right)(x^\ast)
\end{bmatrix}
\end{aligned}
$$
* Then we also have the "stochastic second derivitives" but because all the cross terms are 0, this is just a singleton $`h_{\sigma \sigma}^\ast (x_i)= \frac{\partial^2 h_i}{\partial \sigma^2}(x^\ast)`$. 


We have completely finished describing how to approximate the second order evolution equations. 

The policy function approximations are 

$$
g_i(x) \approx x_i^\ast+g_x^*(x_i)\widehat{x}+\frac{1}{2}\left[ \widehat{x}^Tg_{xx}^\ast(x_i)\widehat{x}+\sigma^2g^\ast_{\sigma\sigma}(x_i) \right]
$$

where are $`g_x^*(x_i),\ g_x^*(x_i), \ g_{\sigma\sigma}^*(x_i)`$ are defined the same way as before but now everything is length $`k`$ 

Recall our earlier structure: $`\textbf{x}^h=(x,h(x))`$:, $`\textbf{y}^g =(g(x),g(h(x)))`$, and 

$$
\begin{aligned}
f(\textbf{x}^h,\textbf{y}^g)&=\begin{bmatrix}
f_1(\textbf{x}^h,\textbf{y}^g) \\
\vdots \\
f_m(\textbf{x}^h,\textbf{y}^g)
\end{bmatrix}
\end{aligned}
$$

Let's say the first $`M`$ equations of $`f`$ correspond to Euler equations. To find the EE erros, we simply plug in our approximation to the first $`M`$ equations. However, this becomes complicated in practice by the fact that we may have multiple state variables. For each state variable, chances are it needs its own grid, unless it's entirely pinned down by other state variables. We have been general up until now, but we can introduce the elastic interest rate example. We can use $`r_{t-1}=r^\ast+\psi[\exp(d_{t-1}-d^\ast)-1]`$ to not have to construct an interest rate grid. Our EE are 

$$
\begin{aligned}
\begin{bmatrix}
 u_c(t)=\beta(1+r_t)\mathbb{E}_t[u_c(t+1)]\\
    -u_h(t)=u_c(t)A_tF_h(t) \\
 u_c(t)\left[1+\Phi'(k_{t+1}-k_t)\right]=\beta \mathbb{E}_t\left[u_c(t+1)\left\{A_{t+1}F_k(t+1)+1-\delta +\Phi'(k_{t+2}-k_{t-1})\right\}\right]
\end{bmatrix}
\end{aligned}
$$

Define the euler equation errors as follows 

$$
\begin{aligned}
ee(\textbf{x},\textbf{y})=\begin{bmatrix}
 -u_c(t)+\beta(1+r_t)\mathbb{E}_t[u_c(t+1)]\\
    u_h(t)+u_c(t)A_tF_h(t) \\
 -u_c(t)\left[1+\Phi'(k_{t+1}-k_t)\right]+\beta \mathbb{E}_t\left[u_c(t+1)\left\{A_{t+1}F_k(t+1)+1-\delta +\Phi'(k_{t+2}-k_{t-1})\right\}\right]
\end{bmatrix}
\end{aligned}
$$

Therefore the errors for our approximation are $` ee(\textbf{x}^{\tilde{h}},\textbf{y}^{\tilde{g}})`$, where $`\tilde{h},\tilde{g}`$ are our second order approximations. To compute the errors. one thing we could do is take the average squared error over grids of capital, technology, and hours. Say that we have $`N `$ points in each grid. Then we could plot the average squared  EE errors vs. capital using the following 

$$
\begin{gather*}
\frac{1}{N^3}\sum_{i=1}^N \sum_{j=1}^N \sum_{l=1}^N ee(\textbf{x}_{ijl}^{\tilde{h}},\textbf{y}_{ijl}^{\tilde{g}})^2 \\
\text{where } \textbf{x}_{ijl}^{\tilde{h}}=(x_{ijl},\tilde{h}(x_{ijl}))=(\textbf{(}k_i,A_j,d_l\textbf{)},\tilde{h}[\textbf{(}k_i,A_j,d_l\textbf{)}]) \\
\text{and } \textbf{y}_{ijl}^{\tilde{g}} =(\tilde{g}(x_{ijl}),\tilde{g}[\tilde{h}(x_{ijl})])=(\tilde{g}\textbf{(}k_i,A_j,d_l\textbf{)},\tilde{g}[\tilde{h}(\textbf{(}k_i,A_j,d_l\textbf{)})])
\end{gather*}
$$

[^1]: I leave $`x_i`$ as an argument since typically our variables have economic meaning which is more useful/intuitive than making row numbers more of a prominent object in interest. To further justify $`x_i`$ as an argument, we consider the following formulation where $`x_i`$ is used to select a row 
