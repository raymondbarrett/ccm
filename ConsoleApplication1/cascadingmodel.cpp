#include <iostream>
#include <vector>
#include <armadillo>
#include <random>
#include "ConsoleApplication1h.h"
#include <experimental/filesystem>
#include <map>



arma::Col<int> step(arma::Col<int>);
void create_res_hall(int[], double[], int, int);
void onetime_load();
arma::Col<int> step(arma::Col<int> currentInfected);
void create_res_hall(int community_sizes[], double community_spread_rates[], int numcommunity, int numMultiSpreaders, std::string path);
arma::Mat<float> create_communities(arma::Mat<float> mat, int numcommunities, int avg_communitysize, float communitysize_variation, float communityweight, float communityweight_variation);
arma::Mat<float> concat_diag_mat(std::experimental::filesystem::path path);
arma::Mat<float> remove_nodes(arma::Mat<float> initial_mat, arma::Col<int> days_to_removal);
void create_dir_graph(int all_community_sizes[][44], double all_community_spread_rates[][44], int numcommunities_per_community[], int numMultiSpreaders_per_community[], std::string res_hall_names[44], std::string path);

int main()
{
	std::cout << "Covid-19 Cascade Simulation\n";

	std::map<std::string, int[]> gatech_apartments = { 
		{"THB",{2,10,60}},
		{"THD", {2,6,20,60}},
		{"THE", {2,6,18,54}},
		{"THF", {2,8,32,96}},
		{"THG", {2,8,24}},
		{"CSN", {4,32,160}},
		{"CSS", {4,48,192}},
		{"CRE", {4,72,360}},
		{"ESE", {4,16,48,192}},
		{"ESS", {4,16,48,192}},
		{"ESW", {4,20,60,240}},
		{"GLC", {4,36,72,360}},
		{"MLD", {4,36,72,288}},
		{"NSL", {4,28,112,448}},
		{"NAE", {4,24,48,624}},
		{"NAN", {4,64,576}},
		{"NAS", {4,36,72,576}},
		{"NAW", {5,45,180}},
		{"ZBR", {4,28,56,224}}
	};
	std::map<std::string, int[]> gatech_suites = 
	std::map<std::string, int[]> gatech_traditional = 
	std::map<std::string, int[]]> greek_housing = 
	

}

arma::Col<int> step(arma::Col<int> currentInfected, arma::Col<int> days_to_removal, int max_days_to_removal)
{
	arma::Mat<float> dir_graph;
	dir_graph.load("processed-mat.txt", arma::arma_binary);
	dir_graph = remove_nodes(dir_graph, days_to_removal);
	currentInfected = arma::conv_to< arma::Col<int> >::from((dir_graph * currentInfected).eval()).row(0);
	for (int i = 0; i < currentInfected.size(); i++)
	{
		if (currentInfected[i] != 0 && currentInfected[i] != 1)
		{
			arma::arma_rng::set_seed_random();
			if (arma::randu<float>() < currentInfected[i]) {
				arma::arma_rng::set_seed_random();
				currentInfected[i] = 1;
				days_to_removal[i] = (int)(arma::randu<float>() * max_days_to_removal);
			}
			
		}
	}
	return currentInfected;
}





/* Subtracts one day off each person who is infected until they have 0 days until removal at which point they are removed from the model
   initial_mat is the initial undirected graph, days_to_removal is a column vector that will be modified
*/
arma::Mat<float> remove_nodes(arma::Mat<float> initial_mat, arma::Col<int> days_to_removal) 
{
	arma::Mat<int> transform (initial_mat.n_cols, initial_mat.n_cols, arma::fill::eye);
	for (int i = 0; i < initial_mat.n_cols; i++) 
	{
		if (days_to_removal[i] > 1) 
		{
			days_to_removal[i] --;
		}
		else if (days_to_removal[i] == 1)
		{
			days_to_removal[i] = 0;
			transform(i, i) = 0;
		}
	}
	return transform * initial_mat * transform;
}


void create_dir_graph(int all_community_sizes[][44], double all_community_spread_rates[][44], int numcommunities_per_community[], int numMultiSpreaders_per_community[], std::string res_hall_names[44], std::string path)
{

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
