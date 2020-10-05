#include <iostream>
#include <vector>
#include <armadillo>
#include <random>
#include "ConsoleApplication1h.h"


arma::Col<int> step(arma::Col<int> currentInfected);
void create_res_hall(int community_sizes[], double community_spread_rates[], int numcommunity, int numMultiSpreaders);
void onetime_load();

int main()
{

	std::cout << "Covid-19 Cascade Simulation";
	int community_sizes[3] = {16,4,2};
	double community_spread_rates[3] = { 0.00032, 0.01, 0.9 };
	create_res_hall(community_sizes, community_spread_rates, 3, 0);
}

arma::Col<int> step(arma::Col<int> currentInfected)
{
	arma::Mat<float> dir_graph;
	dir_graph.load("processed-mat.txt", arma::arma_binary);
	currentInfected = arma::conv_to< arma::Col<int> >::from((dir_graph * currentInfected).eval()).row(0);
	for (int i = 0; i < currentInfected.size(); i++)
	{
		if (currentInfected.at(i) != 0 && currentInfected.at(i) != 1)
		{
			currentInfected.at(i) = arma::randu<int>() < currentInfected.at(i) ? 1 : 0;
		}
	}
	return currentInfected;
}

void create_res_hall(int community_sizes[], double community_spread_rates[], int numcommunity, int numMultiSpreaders)
{
	arma::Mat<float> dir_graph(community_sizes[0] + numMultiSpreaders, community_sizes[0] + numMultiSpreaders, arma::fill::zeros);
	for (int i = 0; i < numcommunity; i++)
	{
		for (int j = 0; j < community_sizes[0]/community_sizes[i]; j++)
		{
			for (int k = j * community_sizes[i]; k < (j+1) * community_sizes[i] - 1; k++)
			{
				for (int p = k + 1; p < (j + 1) * community_sizes[i] - 0; p++)
				{
					dir_graph(k, p) = dir_graph(p, k) = community_spread_rates[i];
				}
			}
		}
	}
	//To implement for multi-spreader code later
	/*for (int q = 0; q < numSuperSpreaders; q++)
	{

	}
	*/
	dir_graph.save("residence_hall_1.txt", arma::raw_ascii);
}

void onetime_load()
{
	arma::Mat<float> mat;
	mat.load("file-name.txt", arma::auto_detect);
	mat.save("processed-mat.txt", arma::arma_binary);
}
