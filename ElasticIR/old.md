
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

Schmitt-Grohé and Uribe made a contribution to the literature by integrating Matlab's symbolic toolbox with this standard economic paradigm. They have a short but well cited [paper](https://faculty.wcas.northwestern.edu/lchrist/papers/perturbation.pdf) that is a good companion for the following discussion, in particular section 4 which goes through an example with 1 control and 2 states. Here, I present the general case and break down each step in their derivation. The basic roadmap is that we will formulate Taylor-type approximations to policy functions and state evolution equations and then measure the quality of these approximations looking at how big the errors in the Euler equations are. 

We denote $` h(x) `$ as the state evolution equations and $`g(x)`$ as the policy functions. To do a second order approximation of the evolution of the $`i\text{-th}`$  state variable, we do the following 

$$
h_i(x) \approx x_i^\ast+h_x^*(x_i)\widehat{x}+\frac{1}{2}\left[ \widehat{x}^Th_{xx}^\ast(x_i)\widehat{x}+\sigma^2h^\ast_{\sigma\sigma}(x_i) \right]
$$

Notation: we have the usual $`x_i^* `$ denoting the steady state of $`x_i`$, $`\widehat{x}_i`$ denoting "linearization" by $`x_i-x_i^* `$, and $` \sigma^2 `$ denoting the variance of the stochastic component of the system. The non-standard notation appear in the coefficient components. These correspond to the coefficients necessary to approximate the solution to the system of equations at the steady state. We will construct the form and notivate the styling of these parts step by step. To begin, we will break down some of the theoretical motication for this procedure. 

Let $`\textbf{x}^h=(x,h(x))`$ and $`\textbf{y}^g =(g(x),g(h(x)))`$. We can exploit the fact that the $`h`$ and $`g`$ must satisfy the following 

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

* To be explicit, $`h`$ is a vector of state evolution equations where $`i\text{-th}`$ element corresponds to $`x_i' = h_i(x)`$, the evolution of the $`i\text{-th}`$ state variable, or

  

$$
\begin{aligned}
h(x)&=\begin{bmatrix}
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

  

* From our original second order equation[^1], $`h_x^*(x_i)`$ is the $`i\text{-th}`$ row of  $`h_x`$ evaluated at the steady state (thus the star).  As seen above, this is all possible (first) partial derivatives of $`h_i`$.

$$
\begin{aligned}
h_x^\ast(x_i)=[h_x^\ast]_{x_i}=\begin{bmatrix}
\frac{\partial h_i}{\partial x_1}(x^\ast) & \cdots & \frac{\partial h_i}{\partial x_n}(x^\ast)
\end{bmatrix}
\end{aligned}
$$

* For the second derivatives, we technically don't have a matrix but instead an array of second derivatives: a collection of $`n`$ matrices of size $`n \times n`$. The first matrix is the partial derivative of the $`h_x`$ matrix w.r.t $`x_1`$, the second matrix are partials w.r.t $`x_2`$ and so on. We are interested in the partial derivatives of $`\nabla h_i`$, and these will just correspond to the $`i\text{-th}`$ row in each matrix of the array. Casting this in linear algebra format  (selecting these rows and then using a row combination; the array we are selecting from is called hxx in Matlab's symbolic toolbox)

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
 \left(\frac{\partial}{\partial x_n}\nabla h_i^T\right)(x^\ast)
\end{bmatrix}
\end{aligned}
$$
* Then we also have the "stochastic second derivitives" but because all the cross terms are 0, this is just a singleton $`h_{\sigma \sigma}^\ast (x_i)= \frac{\partial^2 h_i}{\partial \sigma^2}(x^\ast)`$. 


We have completely finished describing how to approximate the second order evolution equations. 

The policy function approximations are 

$$
g_i(x) \approx x_i^\ast+g_x^*(x_i)\widehat{x}+\frac{1}{2}\left[ \widehat{x}^Tg_{xx}^\ast(x_i)\widehat{x}+\sigma^2g^\ast_{\sigma\sigma}(x_i) \right]
$$

where are $`g_x^*(x_i),\ g_x^*(x_i), \ g_{\sigma\sigma}^*(x_i)`$ are defined the same way as before but now everything is length $`k`$ (still $`n`$ columns). 

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

Let's say the first $`M`$ equations of $`f`$ correspond to Euler equations. To find the EE erros, we simply plug in our approximation to the first $`M`$ equations. We will now use the elastic interest rate example to explicitly illustrate. The EE are 

$$
\begin{aligned}
EE(\textbf{x},\textbf{y})=\begin{bmatrix}
 u_c(t)=\beta(1+r_t)\mathbb{E}\_t[u_c(t+1)]\\
    -u_h(t)=u_c(t)A_tF_h(t) \\
 u_c(t)\left[1+\Phi'(k_{t+1}-k_t)\right]=\beta \mathbb{E}\_t\left[u_c(t+1)\left\\{A_{t+1}F_k(t+1)+1-\delta +\Phi'(k_{t+2}-k_{t-1})\right\\}\right]
\end{bmatrix}
\end{aligned}
$$

Define the Euler equation errors as follows 

$$
\begin{aligned}
ee(\textbf{x},\textbf{y})=\begin{bmatrix}
 -u_c(t)+\beta(1+r_t)\mathbb{E}\_t[u_c(t+1)]\\
    u_h(t)+u_c(t)A_tF_h(t) \\
 -u_c(t)\left[1+\Phi'(k_{t+1}-k_t)\right]+\beta \mathbb{E}\_t\left[u_c(t+1)\left\\{A_{t+1}F_k(t+1)+1-\delta +\Phi'(k_{t+2}-k_{t-1})\right\\}\right]
\end{bmatrix}
\end{aligned}
$$

Therefore the errors for our approximation are $` ee(\textbf{x}^{\tilde{h}},\textbf{y}^{\tilde{g}})`$, where $`\tilde{h},\tilde{g}`$ are our second order approximations.

However, in this example, we have an expectation term. This means we can't simply use a grid to compute the EE errors because we need to take an expectation, which corresponds to the integration of a function. To handle this, we first do a simulation to find the distribution of the economy. This simulation will deliver a large matrix of simulated $`x`$ and $`y`$ points. If we column combine the matrices of  $`x`$ and $`y`$ values, we can think of each row as the economy at a certain point in time $`t`$. For each point in time, we can compute an expectation conditional on today's realization of the stochastic state variable, which also will deliver the EE error for a given  $`t`$ using the equations above. 

To actually calculate these expectations, we use Gauss-Hermite quadrature since the stochastic variable of interest is log-normally distributed. The conditional expectation of $`Z_{t+1}=\ln(A_{t+1})`$ is $`\rho Z_t`$ and the conditional variance is $`\sigma^2`$. For a given order $`n`$, we get $`n`$ nodes $`z_i`$,  use the change of variables formula $`\sigma z_i\sqrt{2}+\rho Z_t`$, and then take the exponential since we have a log-normal. This gives us a quasi-grid of technology nodes to integrate over, with the integration weights $`\omega_i`$ also being determined by the order $`n`$ Gauss-Hermite procedure.   

Putting this all together, to compute the (squared) EE errors at a given $`t`$, with technology nodes $`\tilde{A}_i`$ and a vector ordering $`x=(d,r,k,A)`$, we can do the following, 

$$
\begin{gather*}
ee(\textbf{x}\_{t,}\textbf{y}\_{t})= \left[\sum_{i=1}^n \omega_i\cdot ee(\textbf{x}\_{t}^{\tilde{h}^{A_i}},\textbf{y}\_{t}^{\tilde{g}^{A_i}})\right]^2 \\
\text{where } \tilde{h}^{A_i} = \begin{bmatrix}
 \tilde{h}_1(x)\\
    \tilde{h}_2(x) \\
    \tilde{h}_3(x) \\
 \tilde{A}\_i
\end{bmatrix} \text{ and } \tilde{g}^{A_i}= \tilde{g}(\tilde{h}^{A_i}(x))
\end{gather*}
$$

We can then draw a histogram of EE errors. For the EE which contain stochastic terms, it seems like the errors are non-negligible, with about half above .01. We can normalize by the magnitude of the EE (I use average of LHS and RHS) to get a better sense of how serious they are. We now see an average (squared) error a bit above 1%, which still isn't great. This result is robust to the order of the quadrature. $`n=3`$ does just as well as $`n=30`$, which does makes some sense because a third order quadrature implies an approximate fifth order polynomial, which is a level of non-linearity not to be sneezed at. The non-stochastic EE suggests that at least in most of the code, there aren't any mistakes, but determining tomorrow's states and controls, which depends on the stochastic element could have typos. 


[^1]: I leave $`x_i`$ as an argument since typically our variables have economic meaning which is more useful/intuitive than making row numbers more of a prominent object in interest. To further justify $`x_i`$ as an argument, we consider the following formulation where $`x_i`$ is used to select a row 
