# Implicitly Restarted Shift-and-Invert Arnoldi Method
This projects contains a MATLAB implementation of the one-sided and two-sided implicitly restarted shift-and-invert Arnoldi methods, applicable to linear matrix pencils of the form $(A,B) \in \mathbb{C}^{m \times n}$. The functions `onesided_iram.m` and `twosided_iram.m` are implementations of the one-sided and two-sided methods to a rectangular and square problem respectively. The scripts `main.m`, `main2.m` and `test_filters` contain examples of how to employ these functions in a practical setting.

## One-sided Method

**Syntax**:

[Q,H,Qhist,Hhist] = onesided_iram(Op,b,min_dim,max_dim,rest,filt_method,pencil_size);

**Input arguments:**
* `Op`: Operator defining the matrix vector product $f(v) = (A- \sigma B)^{-1}Bv$
* `b`: Starting vector for the Krylov subspace $\in \mathbb{C}^{n}$
* `min_dim`: defines the number of eigenvalues that are approximated
* `max_dim`: defines the dimension threshold for when a restart occurs
* `rest`: defines the number of restarts the algorithm performs
* `filt_method`: defines how spurious eigenvalues are filtered
  * `0`: Keep largest Ritz values
  * `1`: Keep Ritz values corresponding to vectors that have the smallest residual
* `pencil_size = [m, n]`: defines the size of the base pencil $(A,B) \in \mathbb{C}^{m \times n}$

**Output arguments:**
* `Q,H`: Arnoldi factorization matrices
  * if no breakdown occurs then, $Q(:,end-1)(A- \sigma B)^{-1}B \approx Q H$
  * if breakdown occurs then, $Q(A- \sigma B)^{-1}B \approx Q H$
* `Qhist`: contains the values of Q at each restart
* `Hhist`: contains the values of H at each restart

## Two-sided Method

**Syntax**:

[V,W,H,K] = twosided_iram(A,B,pencil_size,sigma,v1,w1,min,max,restarts,filt_method)

**Input arguments:**
* `A`: First pencil matrix
* `B`: Second pencil matirx
* `pencil_size = [m, n]`: defines the size of the base pencil $(A,B) \in \mathbb{C}^{m \times n}$
* `sigma`: defines the shift towards convergence is prioritized $(\sigma \in \mathbb{C})$
* `v1`: Starting vector for the right Krylov subspace $\in \mathbb{C}^{n}$
* `w1`: Starting vector for the left Krylov subspace $\in \mathbb{C}^{m}$
* `min`: defines the number of eigenvalues that are approximated
* `max`: defines the dimension threshold for when a restart occurs
* `restarts`: defines the number of restarts the algorithm performs
* `filt_method`: defines how spurious eigenvalues are filtered
  * `1`: Keep largest Ritz values
  * `2`: Keep largest Ritz values + Infinite polynomial shift
  * `3`: Keep Ritz values corresponding to vectors that have the smallest residual
  * `4`: Keep Ritz values corresponding to vectors that have the smallest residual + Infinite polynomial shift

**Output arguments:**
* `V,H`: Arnoldi factorization matrices for the Right Krylov Subspace
  * if no breakdown occurs then, $V(:,end-1)(A- \sigma B)^{-1}B \approx V H$
  * if breakdown occurs then, $V(A- \sigma B)^{-1}B \approx V H$
* `W,K`: Arnoldi factorization matrices for the Left Krylov Subspace
  * if no breakdown occurs then, $W(:,end-1)(A- \sigma B)^{-*}B^* \approx W K$
  * if breakdown occurs then, $W(A- \sigma B)^{-*}B^* \approx W K$

