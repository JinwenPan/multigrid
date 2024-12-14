#include "iostream"
#include "memory"
#include "cmath"
#include "fstream"
#include "Grid.h"
#define PI 3.14159265358979

void Grid::set_BC_for_u_Nm(){		// Sets boundary condition for the bonous case (Neumann on vertical edges)
	for (int i=0; i<n; ++i){
		for (int j=0; j<n; ++j){
			u(i,0)=cos(PI*h*i);		// Dirichlet on bottom boundary
			u(i,n-1)=-cos(PI*h*i);	// Dirichlet on top boundary
		}
	}
}

void Grid::set_BC_for_u(){		// Sets boundary condition for the base case (all Dirichlet)
	for (int i=0; i<n; ++i){
		for (int j=0; j<n; ++j){
			u(0,j)=cos(PI*h*j);		// left boundary
			u(n-1,j)=-cos(PI*h*j); 	// right boundary
			u(i,0)=cos(PI*h*i);		// bottom boundary
			u(i,n-1)=-cos(PI*h*i);	// top boundary
		}
	}
}

void Grid::initialise_b(){		// Initialise the b array
	for (int i=0; i<n; ++i){
		for (int j=0; j<n; ++j){
			b(i,j)= 2*PI*PI*cos(PI*h*i)*cos(PI*h*j);
		}
	}
}

void Grid::set_solution(){		// Sets b to the analyticial solution (for finding error)
	for (int i=0; i<n; ++i){
		for (int j=0; j<n; ++j){
			b(i,j)=cos(PI*h*i)*cos(PI*h*j);
		}
	}
}

void Grid::compute_r(){		// Computes the residual from u and b
	double C = 1/(h*h);

	for (int i=1; i<n-1; ++i){ 		// Iterates over iterior nodes (error over boundary is always zero)
			for (int j=1; j<n-1; ++j){
				r(i,j) =  C * (4*u(i,j) - u(i-1,j) - u(i+1,j) - u(i,j-1) - u(i,j+1)) - b(i,j);
			}
	}
}

void Grid::relax(int iter){		// Gauss Seidel scheme to solve for u. 'iter' argument is the number of sweeps
	double C = (h*h)/4;

	for (int it=0; it<iter; ++it){
		for (int i=1; i<n-1; ++i){		// Iterates over interior nodes (u over all boundaries is set by the Dirichlet boundary condtion)
			for (int j=1; j<n-1; ++j){
				u(i,j) = C * b(i,j) + 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1));
			}
		}
		
	}
}

void Grid::relax_Nm(int iter){	// Gauss Seidel solver for the Neumann case
	double C = (h*h)/4;

	for (int it=0; it<iter; ++it){
		for (int i=1; i<n-1; ++i){		// Iterates over interior nodes but left and right Neumann boundaries are treated seperately (top and bottom edges are untouched)
			for (int j=1; j<n-1; ++j){

				u(0,j) = C * b(0,j) + 0.25*(2*u(1,j) + u(0,j-1) + u(0,j+1));	// u on the left boundary
				u(n-1,j) = C * b(n-1,j) + 0.25*(2*u(n-2,j) + u(n-1,j-1) + u(n-1,j+1));	// u on the right boundary
				u(i,j) = C * b(i,j) + 0.25*(u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1));	// u on the rest of interior nodes
			}
		}
	}
}

void Grid::set_u_to_zero(){		// Sets u to zero (when moving down to coarser grid in the v-cycle, 
								// u is set to zero because it is actually the error that needs to be found on the coarse grid)
	for (int i=1; i<n-1; ++i){	// Loops over iterior nodes because error on the boundary is always zero
		for (int j=1; j<n-1; ++j){
			u(i,j)= 0.0;
		}
	}
}

void Grid::restrict_r(Grid &omega2){	// Restricts residual from current grid to the coarser grid (omega2) which is passed in as an argument
	int n2 = omega2.n;	// n of the coarse grid
	auto b2 = omega2.b;	// b of the coarse grid
	for (int i=1; i<n2-1 ; ++i){	// The residual from the current fine grid is restricted on the b array for the coarse grid (omega2)
		for (int j=1; j<n2-1 ; ++j){
			b2(i,j) = 0.125 *(r((2*i-1), 2*j-1) + r((2*i-1), 2*j+1) + r((2*i+1), +2*j-1) + r((2*i+1), 2*j+1))
										+ 0.25 * ( r((2*i), 2*j-1) + r((2*i), 2*j+1) + r((2*i-1), 2*j) + r((2*i+1), 2*j))
										+ 0.5 * r((2*i), 2*j);
		}
	}
}

void Grid::prolong_u(Grid &omega2){ // Prolongs u to the finer grid (omega2) which is passed in as an argument
	auto e2 = omega2.e;	// e of the fine grid
	for (int i=0; i<n-1 ; ++i){	// u from the current coarse grid (i.e. the error) is prolonged to the error array for the fine grid (omega2)
		for (int j=0; j<n-1 ; ++j){
			e2((2*i), 2*j) =   0.5 * u(i, j);
			e2((2*i+1), 2*j) =  0.25*(u(i, j) + u((i+1), j));
			e2((2*i), 2*j+1) =  0.25*(u(i, j) + u(i,j+1));
			e2((2*i+1), 2*j+1) =  0.125*(u(i, j) + u((i+1), j) + u(i, j+1) + u((i+1), j+1));
		}
	}
}

void Grid::correct_u(){ // The error prolonged from the coarse grid to the current fine grid (saved in e) is used to correct u on the current grid
	for (int i=1; i<n-1 ; ++i){
		for (int j=1; j<n-1 ; ++j){
			u(i, j) -= e(i, j);
		}
	}
}

double Grid::L2_norm(const Data arr) const{ // Const function that takes an array of type Data and returns the discrete L2 norm as a double
	double square_sum=0.0;
	for (int i=0; i<n; ++i){
		for (int j=0; j<n; ++j){
			square_sum += arr(i,j)*arr(i,j);
		}
	}
	return sqrt(square_sum/grid_size);
}

void Grid::compute_e(Grid &analytic_sol){ // Computes the error at each node w.r.t the analytical solution
	for (int i=0; i<n; ++i){
		for (int j=0; j<n; ++j){
			e(i,j) = u(i,j) - analytic_sol.b(i,j);
		}
	}
}

void Grid::write_output() const{	// Function to write soluton u to a .txt file for gnuplot
	std::ofstream GS_output;
	GS_output.open("solution.txt");
	std::cout<< "Wrting output " << std::flush;
	GS_output << "#X\tY\tu" << std::endl;
	for (int i=0; i<n; ++i){
			for (int j=0; j<n; ++j){
				GS_output << i*h << "\t" << j*h << "\t" << u(i,j) << std::endl;
			}
	}
	std::cout << "done." << std::endl;
}
