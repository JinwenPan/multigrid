#ifndef GRID_H
#define GRID_H
#include "Data.h"

struct Grid{	// This structure holds all the data associated with each level and perfoms functions on them (all members are public)

	Grid()= default;	// Default constructor
	Grid(int &N) : n(N), u(N), b(N), r(N), e(N){};	// N is the number of grid points in either x or y direction. 
		// N is the only argument passed in when instantiating a Grid object. All the data is constructed from N
	__attribute__ ((noinline)) ~Grid() = default;	// Default destructor forced to not be inlined by compiler
	int n;	// n is the number of grid points in either x or y direction
	int grid_size = n*n;	// Total number of grid points
	double h = 1.0/(n-1);	// Vertical or horizontal spacing between nodes
	Data u, b, r, e;		// Data on the current grid: u(solution), b (right hand side), r(residual) and e(error)

	// Description of function members below are provided in the implementation code in Grid.cpp
	void set_BC_for_u();
	void set_BC_for_u_Nm();
	void initialise_b();
	void set_solution();
	void relax(int);
	void relax_Nm(int);
	void set_u_to_zero();
	void compute_r();
	void compute_e(Grid&);
	void restrict_r(Grid&);
	void prolong_u(Grid&);
	void correct_u();
	double L2_norm(const Data) const;
	void write_output() const;
	
};

#endif