#pragma once
#include <armadillo>

class cascadingmodel
{
public:
	arma::Col<int> step(arma::Col<int>);
	void create_res_hall(int[], double[], int, int);
	void onetime_load();
};