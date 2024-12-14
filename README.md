# Multigrid Method

This repository provides a code implementation to solve the 2D Laplace equation in Cartesian coordinates using the finite difference method. The linear system generated by discretizing the equation is solved using the V-cycle multigrid method. Both Dirichlet and Neumann boundary conditions are supported. L2 norm is used to compute the residual, error, and convergence factor.

## Requirements

Ensure you have the necessary environment to compile and run the program. The code uses a standard `Makefile` to build the executable.

### Compilation

To compile the program, save all files in one directory and run the following command:

```bash
make
```

This will compile the program and generate the executable named `mgsolve`.

## Usage

The program expects two arguments (for Dirichlet boundary condition on the edges of [0, 1]x[0, 1]):
1. The number of levels for the multigrid method.
2. The number of V-cycle iterations.

Run the program as follows:

```bash
./mgsolve 7 20
```

where:
- `7` is the number of multigrid levels.
- `20` is the number of V-cycle iterations.

### Neumann on Vertical Edges

To activate the Neumann boundary condition on vertical edges, provide a third argument. For example:

```bash
./mgsolve 7 130 n
```

where `n` indicates the activation of the bonus task.

## Output

The program will solve the Laplace equation and output the results, including the residual, error, and convergence factor in the terminal and numerical solution in `solution.txt`.
