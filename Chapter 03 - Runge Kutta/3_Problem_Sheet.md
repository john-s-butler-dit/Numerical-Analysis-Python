<!-- #region -->
# Problem Sheet 3 - Runge Kutta
1. a) Apply the Midpoint Method to approximate the solution of the given initial value problems using the indicated number of time steps. Compare the approximate solution with the given exact solution:
$$ y'=t-y, \ \ (0\leq t \leq 4),$$
with the initial condition $y(0)=1,$
Let $N=4$,  with the exact solution
$$y(t)=2e^{-t}+t-1.$$


1. b) Apply the Midpoint Method to approximate the solution of the given initial value problems using the indicated number of time steps. Compare the approximate solution with the given exact solution:
$$y'=y-t, \ \ (0\leq t \leq 2),$$
with the initial condition $y(0)=2,$
Let $N=4$, with the exact solution
$$y(t)=e^{t}+t+1.$$



2. a) Apply the 4th Order Runge Kutta Method to approximate the solution of the given initial value problems using the indicated number of time steps. Compare the approximate solution with the given exact solution
$$y'=t-y, \ \ (0\leq t \leq 4),$$
with the initial condition $y(0)=1,$
Let $N=4$, with the exact solution
$$y(t)=2e^{-t}+t-1.$$


2. b) Apply the 4th Order Runge Kutta Method to approximate the solution of the given initial value problems using the indicated number of time steps. Compare the approximate solution with the given exact solution
$$y'=y-t, \ \ (0\leq t \leq 2)$$
with the initial condition $y(0)=2,$
$N=4$, with the exact solution
$$y(t)=e^{t}+t+1.$$



3. Derive the difference equation for the Midpoint Runge Kutta method
$$ w_{n+1}=w_n+k_2,$$
$$k_1=hf(t_n,w_n),$$
$$k_2=hf(t_n+\frac{1}{2}h,w_n+\frac{1}{2}k_1),$$
for solving the ordinary differential equation
$$ \frac{dy}{dt}=f(t,y), $$
$$ y(t_0)=y_0, $$
by using a formula of the form
$$ w_{n+1}=w_n+ak_1+bk_2, $$
where $k_1$ is defined as above,
$$k_2=hf(t_n+\alpha h,w_n+\beta k_1),$$
and $a$, $b$, $\alpha$ and $\beta$ are constants are determined. Prove that $a+b=1$ and $b\alpha=b\beta=\frac{1}{2}$ and choose appropriate values to give the Midpoint Runge Kutta method.


<!-- #endregion -->
```python

```
