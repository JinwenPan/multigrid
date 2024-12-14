#ifndef DATA_H
#define DATA_H
#include "memory"

struct Data{	// This struct holds an array of elements on dynamic memory and overloads operator () for easy accessibility (all members are public)

	Data()= default;
	Data(int N) : n(N), arr(new double[n*n]){} // Constructs n and an arr from the N variable passed to a Data object

	inline double& operator ()(int i, int j) {return arr[i*n +j];} 		// Returns Lvalue (overloads () and returns a reference for modifying the (i-th, j-th) element)
	inline double operator ()(int i, int j) const {return arr[i*n +j];} 	// Returns Rvalue (overloads () as a read-only member function for getting the (i-th, j-th) element)

	int n;	// Number of elements in either x or y direction
	std::shared_ptr<double[]> arr;	// arr declared as a shared_ptr so the memory is appropriately freed without the need for an explicit destructor
};

#endif