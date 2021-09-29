<!-- #region -->
# Problem Sheet 4 - Multistep Methods
1. a) Apply the 3-step Adams-Bashforth to approximate the solution of the given initial value problems using the indicated number of time steps. Compare the approximate solution with the given exact solution:
$$ y'=t-y, \ \ (0\leq t \leq 4),$$
with the initial condition $y(0)=1,$
Let $N=4$,  with the exact solution
$$y(t)=2e^{-t}+t-1.$$


1. b) Apply the 3-step Adams-Bashforth to approximate the solution of the given initial value problems using the indicated number of time steps. Compare the approximate solution with the given exact solution:
$$y'=y-t, \ \ (0\leq t \leq 2),$$
with the initial condition $y(0)=2,$
Let $N=4$, with the exact solution
$$y(t)=e^{t}+t+1.$$

2. a) Apply the 2-step Adams-Moulton Method to approximate the solution of the given initial value problems using the indicated number of time steps. Compare the approximate solution with the given exact solution:
$$ y'=t-y, \ \ (0\leq t \leq 4),$$
with the initial condition $y(0)=1,$
Let $N=4$,  with the exact solution
$$y(t)=2e^{-t}+t-1.$$

2. b) Apply the 2-step Adams-Moulton Method to approximate the solution of the given initial value problems using the indicated number of time steps. Compare the approximate solution with the given exact solution:
$$y'=y-t, \ \ (0\leq t \leq 2),$$
with the initial condition $y(0)=2,$
Let $N=4$, with the exact solution
$$y(t)=e^{t}+t+1.$$



3. Derive the difference equation for the 1-step Adams-Bashforth method:
$$ w_{n+1}=w_n+hf(t_{n},w_{n}),$$
with the local truncation error
$$ \tau_{n+1}(h)=\frac{h}{2}y^{2}(\mu_n),$$
where $\mu_n \in (t_{n},t_{n+1})$.

4. Derive the difference equation for the 2-step Adams-Bashforth method:
$$ w_{n+1}=w_n+(\frac{3}{2}hf(t_{n},w_{n})-\frac{1}{2}hf(t_{n-1},w_{n-1})),$$
with the local truncation error
$$ \tau_{n+1}(h)=\frac{5h^2}{12}y^{3}(\mu_n),$$
where $\mu_n \in (t_{n-1},t_{n+1})$.

5.  Derive the difference equation for the 3-step Adams-Bashforth method:
$$ w_{n+1}=w_n+(\frac{23}{12}hf(t_{n},w_{n})-\frac{4}{3}hf(t_{n-1},w_{n-1})+\frac{5}{12}hf(t_{n-2},w_{n-2})),$$
with the local truncation error
$$ \tau_{n+1}(h)=\frac{9h^3}{24}y^{4}(\mu_n),$$
where $\mu_n \in (t_{n-2},t_{n+1})$.

6. Derive the difference equation for the 0-step Adams-Moulton method:
$$ w_{n+1}=w_n+hf(t_{n+1},w_{n+1}),$$
with the local truncation error
$$\tau_{n+1}(h)=-\frac{h}{2}y^{2}(\mu_n),$$
where $\mu_n \in (t_{n-2},t_{n+1})$.

7. Derive the difference equation for the 1-step Adams-Moulton method:
$$ w_{n+1}=w_n+\frac{1}{2}hf(t_{n+1},w_{n+1})+\frac{1}{2}hf(t_{n},w_{n}),$$
with the local truncation error
$$ \tau_{n+1}(h)=-\frac{h^2}{12}y^{3}(\mu_n),$$
where $\mu_n \in (t_{n},t_{n+1})$.

8. Derive the difference equation for the 2-step Adams-Moulton method:
$$ w_{n+1}=w_n+\frac{5}{12}hf(t_{n+1},w_{n+1})+\frac{8}{12}hf(t_{n},w_{n})-\frac{1}{12}hf(t_{n-1},w_{n-1}),$$
with the local truncation error
$$ \tau_{n+1}(h)=-\frac{h^3}{24}y^{4}(\mu_n),$$
where $\mu_n \in (t_{n-1},t_{n+1})$.
9.  Derive the difference equation for the 3-step Adams-Moulton method:\\
$$ w_{n+1}=w_n+\frac{9}{24}hf(t_{n+1},w_{n+1})+\frac{19}{24}hf(t_{n},w_{n})-\frac{5}{24}hf(t_{n-1},w_{n-1})+\frac{1}{24}hf(t_{n-2},w_{n-2}),$$
with the local truncation error
$$ \tau_{n+1}(h)=-\frac{h^4}{720}y^{5}(\mu_n),$$
where $\mu_n \in (t_{n-2},t_{n+1})$.


<!-- #endregion -->
```python

```
