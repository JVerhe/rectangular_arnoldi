# Implicitly Restarted Shift-and-Invert Arnoldi Method
This projects contains a MATLAB implementation of the one-sided and two-sided implicitly restarted shift-and-invert Arnoldi methods. The scripts `main.m` and `main2.m` apply the one-sided and two-sided methods to a rectangular and square problem respectively.

## One-sided Method

The script `main.m` calls the one-sided Arnoldi method to solve a rectangular pencil defined by the parameters in the first lines of the program. The parameters are the following:

* `m, n`: define the size of the base pencil $(A,B) \in \mathbb{C}^{m \times n}$
* `sigma`: defines the shift towards convergence is prioritized $(\sigma \in \mathbb{C})$
* `rest`: defines the number of restarts the algorithm performs
* `min_dim`: defines the number of eigenvalues that are approximated
* `max_dim`: defines the dimension threshold for when a restart occurs
* `filt_method`: defines how spurious eigenvalues are filtered
  * `1`: converge towards largest ritz values
  * `2`: converge towards ritz vectors with smallest residual part
* `pencil_type`: defines the type of rectangular pencil
  * `1`

## Two-sided Method