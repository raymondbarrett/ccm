#include <iostream>
#include <vector>
#include <armadillo>
#include <random>
#include "ConsoleApplication1h.h"
#include <experimental/filesystem>


arma::Col<int> step(arma::Col<int>);
void create_res_hall(int[], double[], int, int);
void onetime_load();
arma::Col<int> step(arma::Col<int> currentInfected);
void create_res_hall(int community_sizes[], double community_spread_rates[], int numcommunity, int numMultiSpreaders, std::string path);
arma::Mat<float> create_communities(arma::Mat<float> mat, int numcommunities, int avg_communitysize, float communitysize_variation, float communityweight, float communityweight_variation);
arma::Mat<float> concat_diag_mat(std::experimental::filesystem::path path);


int main()
{

	std::cout << "Covid-19 Cascade Simulation";
	/*int community_sizes[3] = {16,4,2};
	int community_sizes_2[4] = { 24,12,6,3 };
	double community_spread_rates_2[4] = { 0.001, 0.02, 0.1, 0.8 };
	double community_spread_rates[3] = { 0.00032, 0.01, 0.9 };
	create_res_hall(community_sizes, community_spread_rates, 3, 0, "graphs/residence_hall_1.txt");
	create_res_hall(community_sizes_2, community_spread_rates_2, 4, 0, "graphs/residence_hall_2.txt");
	*/
	concat_diag_mat(std::experimental::filesystem::path("graphs")).save("graphs/residence_hall_concat.txt", arma::raw_ascii);
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




/*  Creates an nxn undirected graph where n is the size of the residence hall. Takes community_sizes, an array where the first entry is
	the largest community and successive entries describe the sizes of smaller subcommunities; in residence halls the entire hall might be the
	largest community, and a certain floor is a subcommunity, with roommates being the smallest subcommunity.
	Community spread rates describes the spread rates each of the subcommunities, community_spread_rates[i] is the chance that if someone in your
	subcommunity described by community_sizes[i], that you get infected. 
	numcommunity is the size of the arrays community_sizes and create_res_hall (because c++ is dumb for some reason)
	numMultiSpreaders not implemented yet, keep at 0. Hopefully will later describe the number of people who interact with all people in a certain residence hall (RA's, dorm heads, etc.)
*/
void create_res_hall(int community_sizes[], double community_spread_rates[], int numcommunity, int numMultiSpreaders, std::string path)
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
	dir_graph.save(path, arma::raw_ascii);
}





/*  numcommunities is the number of disparate communities that the graph has, assume that the size of these communities follows a normal distribution centered at
	avg_communitysize with standard deviation communitysize_variation and the weight of these connections individually being another normal distribution centered at
	communityweight and with standard deviation communityweight_variation.
*/
//use 0.1 for currentcommunityweight
arma::Mat<float> create_communities(arma::Mat<float> mat, int numcommunities, int avg_communitysize, float communitysize_variation, float communityweight, float communityweight_variation) 
{
	std::normal_distribution <float> communitysize_dist(avg_communitysize, communitysize_variation);
	std::normal_distribution <float> communityweight_dist(communityweight, communityweight_variation);
	std::default_random_engine generator;
	int currentCommunitySize;
	int currentCommunityWeight;
	int pop = mat.n_cols;
	for (int i = 0; i < numcommunities; i++)
	{
		currentCommunitySize = int(communitysize_dist(generator));
		currentCommunityWeight = communitysize_dist(generator);
		std::vector<int> randarray(pop);
		for (int j = 0; j < pop; j++) 
		{
			randarray[i] = i;
		}
		std::shuffle(std::begin(randarray), std::end(randarray), generator);


		for (int k = 0; k < currentCommunitySize; k++) 
		{
			for (int m = k + 1; m < currentCommunitySize; m++) 
			{
				mat(randarray[k], randarray[m]) = mat(randarray[m], randarray[k]) = mat(randarray[k], randarray[m]) + currentCommunityWeight;
			}
		}
	}
	return mat;
}





//Concatenates all matrices in file path mat_path into one output matrix loaded into memory
arma::Mat<float> concat_diag_mat(std::experimental::filesystem::path path)
{
	namespace fsys = std::experimental::filesystem;
	const fsys::directory_iterator end{};
	arma::Mat<float> concat_mat;
	arma::Mat<float> cache;
	for (fsys::directory_iterator mat{ path }; mat != end; mat++) 
	{
		cache.load(mat->path().string(), arma::auto_detect);
		int init_dim = concat_mat.n_cols;
		concat_mat.resize(init_dim + cache.n_cols, init_dim + cache.n_cols);
		for (int i = init_dim; i < cache.n_cols; i++)
		{
			for (int j = i; j < cache.n_cols; j++) 
			{
				concat_mat(i, j) = concat_mat(j, i) = cache(i - init_dim, j - init_dim);
			}
		}
		
	}
	return concat_mat;
}



//Converts a .txt matrix of any format and converts it into arma_binary format
void onetime_load()
{
	arma::Mat<float> mat;
	mat.load("../Resource Files/file-name.txt", arma::auto_detect);
	mat.save("../Resource Files/processed-mat.txt", arma::arma_binary);
}
