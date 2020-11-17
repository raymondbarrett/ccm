#include <iostream>
#include <vector>
#include <armadillo>
#include <random>
#include "ConsoleApplication1h.h"
#include <experimental/filesystem>
#include <map>



arma::Col<float> step(arma::Mat<float> dir_graph, arma::Col<int> currentInfected, arma::Col<int> days_to_removal, int max_days_to_removal);
void onetime_load();
void create_res_hall(std::vector<int> community_sizes, std::vector<double> community_spread_rates, int numMultiSpreaders, std::string path);
arma::Mat<float> create_communities(arma::Mat<float> mat, int numcommunities, int avg_communitysize, float communitysize_variation, float communityweight, float communityweight_variation);
arma::Mat<float> concat_diag_mat(std::experimental::filesystem::path path);
int create_dir_graph(std::map<std::string, float> community_info, std::map<std::string, std::vector<int>> apts, std::map<int, std::vector<double>> apts_rates, std::map<std::string, std::vector<int>> suite, std::map<int, std::vector<double>> suite_rates, std::map<std::string, std::vector<int>>trad, std::map<int, std::vector<double>> trad_rates, std::map<std::string, std::vector<int>> greek, std::map<int, std::vector<double>> greek_rates, std::string path);
void step_manager(arma::Col<int> currentInfected, arma::Col<int> days_to_removal, int max_days_to_removal, std::string path, int numDays);


int main()
{
	std::cout << "Covid-19 Cascade Simulation\n";
	//create_res_hall({ 2,10,60 }, { .97, .02, .001 }, 0, "graphs/test.txt");
	std::map<std::string, std::vector<int>> gatech_apartments = { 
		{"THB",	{60,2,10}},
		{"THD", {60,20,6,2}},
		{"THE", {54,18,6,2}},
		{"THF", {96,32,8,2}},
		{"THG", {24,8,2}},
		{"CSN", {160,32,4}},
		{"CSS", {192,48,4}},
		{"CRE", {360,72,4}},
		{"ESE", {192,48,16,4}},
		{"ESS", {192,48,16,4}},
		{"ESW", {240,60,20,4}},
		{"GLC", {360,72,36,4}},
		{"MLD", {288,72,36,4}},
		{"NSL", {448,112,28,4}},
		{"NAE", {624,48,24,4}},
		{"NAN", {576,64,4}},
		{"NAS", {576,72,36,4}},
		{"NAW", {180,45,5}},
		{"ZBR", {224,56,28,4}}
	};
	std::map < int, std::vector<double>> apartment_rates =
	{
		{3, {.00001,.0002,.75}},
		{4, {.00001,.0001,.0004,.75}}
	};
	std::map<std::string, std::vector<int>> gatech_suites = {
		{"WDS", {128,32,16,2}},
		{"GLDSTN", {48,24,12,2}},
		{"HAYGRY", {48,24,12,2}},
		{"HRS", {96,32,16,4,2}},
		{"WDN", {128,32,16,4,2}},
	};
	std::map < int, std::vector<double>> suite_rates =
	{
		{4, {.00001,.0001,.0005,.75}},
		{5, {.00001,.0001,.0005,.75,.95}}
	};
	std::map<std::string, std::vector<int>> gatech_traditional = {
		{"ARM", {114,38,2}},
		{"BRN", {90,30,2}}, 
		{"CAL", {152,38,2}},
		{"CLD", {120,40,20,2}},
		{"FLD", {128,32,2}},
		{"FIT", {128,32,2}},
		{"FLK", {160,40,2}},
		{"FRE", {114,38,2}},
		{"FUL", {60,30,2}},
		{"GLN", {360,72,36,2}},
		{"HAN", {120,30,2}},
		{"HRN", {150,30,2}},
		{"HEF", {114,38,2}},
		{"HOP", {128,32,2}},
		{"HOW", {126,42,2}},
		{"MTH", {152,38,2}},
		{"MON", {114,38,2}},
		{"PER", {120,30,2}},
		{"SMT", {288,72,36,2}},
		{"TOW", {256,64,32,2}}
	};
	std::map < int, std::vector<double>> traditional_rates =
	{
		{3, {.00001,.01,.95}},
		{4, {.00001,.0001,.01,.95}}
	};
	std::map < std::string, std::vector<int>> greek_housing = {
		{"PKA", {32}},
		{"PKP", {32}},
		{"PKTh", {32}},
		{"SPE", {32}},
		{"PGD", {32}},
		{"KS", {32}},
		{"ZBT", {32}},
		{"TC", {32}},
		{"PU", {32}},
		{"PSK", {32}},
		{"TKE", {32}},
		{"ASP", {32}},
		{"DC", {32}},
		{"DU", {32}},
		{"KA", {32}},
		{"PKS", {32}},
		{"TX", {32}},
		{"SC", {32}},
		{"DSP", {32}},
		{"PKTa", {32}},
		{"ATO", {32}},
		{"SAE", {32}},
		{"DTD", {32}},
		{"BTP", {32}},
		{"SN", {32}},
		{"PDT", {32}},
		{"CPh", {32}},
		{"AEP", {32}},
		{"CPs", {32}},
		{"ADC", {32}},
		{"ACO", {32}},
		{"ADP", {32}},
		{"AP", {32}},
		{"PM", {32}},
		{"ZTA", {32}},
		{"AGD", {32}},
		{"AXD", {32}},
		{"KAT", {32}}
	};
	std::map < int, std::vector<double>> greek_housing_rates = 
	{
		{1, {.6}}
	};
	std::map <std::string, float> community_info = 
	{
		{"numcommunities", 500},
		{"avg_communitysize", 8},
		{"communitysize_variation", 3},
		{"avg_communityweight", 0.05},
		{"communityweight_variation", 0.3}
	};
	int num_seeds = 3;
	int max_days_to_removal = 12;

	int mat_size = create_dir_graph(community_info, gatech_apartments, apartment_rates, gatech_suites, suite_rates, gatech_traditional, traditional_rates, greek_housing, greek_housing_rates, "big_mat.txt");
	//int mat_size = 72;
	std::vector<int> shuffle(mat_size);
	for (int i = 0; i < shuffle.size(); i++) 
	{
		shuffle[i] = i;
	}
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine rand_engine(seed);
	std::shuffle(std::begin(shuffle), std::end(shuffle), rand_engine);
	arma::Col<int> currentInfected(mat_size, arma::fill::zeros);
	arma::Col<int> days_to_removal(mat_size, arma::fill::zeros);
	std::uniform_int_distribution <int> removal_time(1, max_days_to_removal);
	for (int j = 0; j < num_seeds; j++) 
	{
		//std::cout << shuffle[j];
		//std::cout << "\n";
		currentInfected[shuffle[j]] = 1;
		days_to_removal[shuffle[j]] = removal_time(rand_engine);
	}
	days_to_removal[0] = 5;
	currentInfected[mat_size] = 1;
	days_to_removal[mat_size] = INT_MAX;
	arma::Col<int> initial_spread;
	std::cout << "\nSimulating Spread...";
	step_manager(currentInfected, days_to_removal, max_days_to_removal, "big_mat.txt", 40);
}

void step_manager(arma::Col<int> currentInfected, arma::Col<int> days_to_removal, int max_days_to_removal, std::string path, int numDays)
{
	arma::Col<int> infected_vec = currentInfected;
	int currentSum = 0;
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine rand_engine(seed);
	std::uniform_int_distribution<int> removaltime(4, max_days_to_removal);
	std::uniform_real_distribution<float> rand_chance(0, 1);
	for (int j = 0; j < infected_vec.n_rows; j++)
	{
		currentSum += infected_vec[j];
	}
	std::cout << "\nNumber of Infected People: ";
	std::cout << currentSum;
	std::string header = "\nCurrently simulating day ";
	arma::Mat<float> dir_graph;
	dir_graph.load(path, arma::auto_detect);
	arma::Row<float> zero_row(dir_graph.n_cols, arma::fill::zeros);
	arma::Col<float> zero_col(dir_graph.n_cols, arma::fill::zeros);
	arma::Col<float> infection_chances(dir_graph.n_cols, arma::fill::zeros);
	//Removes nodes that are no longer infectious
	for (int i = 0; i < numDays; i++) {
		std::cout << header + std::to_string(i + 1);
		
		for (int j = 0; j < days_to_removal.size(); j++) 
		{
			if (days_to_removal[j] > 1) 
			{
				days_to_removal[j] = days_to_removal[j] - 1;
			}
			else if (days_to_removal[j] == 1) 
			{
				days_to_removal[j] = 0;
				infected_vec[j] = 0;
				dir_graph.row(i) = zero_row;
				dir_graph.col(i) = zero_col;
			}
		}
		//std::cout << "\n";
		/*for (int k = 0; k < days_to_removal.size(); k++) 
		{
			std::cout << std::to_string(days_to_removal[k]);
		}
		*/
		//std::cout << "\n";
		infection_chances = step(dir_graph, infected_vec, days_to_removal, max_days_to_removal);
		for (int m = 0; m < infection_chances.size(); m++) 
		{
			if (rand_chance(rand_engine) < infection_chances[m])
			{
				infected_vec[m] = 1;
				days_to_removal[m] = removaltime(rand_engine);
			}
		}

		currentSum = 0;
		for (int j = 0; j < infected_vec.n_rows; j++)
		{
			currentSum += infected_vec[j];
		}
		std::cout << "\nNumber of Infected People: ";
		std::cout << currentSum;
	}
}


arma::Col<float> step(arma::Mat<float> dir_graph, arma::Col<int> currentInfected, arma::Col<int> days_to_removal, int max_days_to_removal)
{
	arma::Col<float> currentInfectedFloat = arma::conv_to < arma::Col<float>>::from(currentInfected);
	arma::Col<float> infection_chances = arma::conv_to<arma::Col<float>>::from( ((dir_graph * currentInfected).eval()).col(0));
	return infection_chances;
}


int create_dir_graph(std::map<std::string, float> community_info, std::map<std::string, std::vector<int>> apts, std::map<int, std::vector<double>> apts_rates, std::map<std::string, std::vector<int>> suite, std::map<int, std::vector<double>> suite_rates, std::map<std::string, std::vector<int>>trad, std::map<int, std::vector<double>> trad_rates, std::map<std::string, std::vector<int>> greek, std::map<int, std::vector<double>> greek_rates, std::string path)
{
	std::string folder = "graphs/";
	std::string ext = ".txt";
	std::string filler = "\nCreating matrix ";
	for (auto const &apt_names : apts) 
	{
		std::cout << filler + apt_names.first;
		create_res_hall(apt_names.second, apts_rates[apt_names.second.size()], 0, folder + apt_names.first + ext);
	}
	for (auto const &suite_names : suite)
	{
		std::cout << filler + suite_names.first;
		create_res_hall(suite_names.second, suite_rates[suite_names.second.size()], 0, folder + suite_names.first + ext);
	}
	for (auto const &trad_names : trad)
	{
		std::cout << filler + trad_names.first;
		create_res_hall(trad_names.second, trad_rates[trad_names.second.size()], 0, folder + trad_names.first + ext);
	}
	for (auto const &greek_names : greek)
	{
		std::cout << filler + greek_names.first;
		create_res_hall(greek_names.second, greek_rates[greek_names.second.size()], 0, folder + greek_names.first + ext);
	}
	
	std::cout << "\nCompressing into one matrix...";
	arma::Mat<float> dir_graph = concat_diag_mat("graphs");
	std::cout << "\nCompressed!";
	std::cout << "\nCreating communities...";
	dir_graph = create_communities(dir_graph, community_info["numcommunities"], community_info["avg_communitysize"], community_info["communitysize_variation"], community_info["avg_communityweight"], community_info["communityweight_variation"]);
	dir_graph.resize(dir_graph.n_cols + 1, dir_graph.n_cols + 1);
	std::normal_distribution <float> ATLtoCampus_weight(0.00005203571, 0.00001);
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();;
	std::default_random_engine generator(seed);
	float rand;
	for (int i = 0; i++; i < dir_graph.n_cols) 
	{
		rand = ATLtoCampus_weight(generator);
		dir_graph(i, dir_graph.n_cols) = dir_graph(dir_graph.n_cols, i) = rand < 0 ? 0 : rand > 1 ? 1 : rand;
	}
	std::cout << "\nCreated!";
	std::cout << "\nSaving...";
	dir_graph.save(path);
	std::cout << "\nSaved!";
	return dir_graph.n_cols;
}


/*  Creates an nxn undirected graph where n is the size of the residence hall. Takes community_sizes, an array where the first entry is
	the largest community and successive entries describe the sizes of smaller subcommunities; in residence halls the entire hall might be the
	largest community, and a certain floor is a subcommunity, with roommates being the smallest subcommunity.
	Community spread rates describes the spread rates each of the subcommunities, community_spread_rates[i] is the chance that if someone in your
	subcommunity described by community_sizes[i], that you get infected. 
	numcommunity is the size of the arrays community_sizes and create_res_hall (because c++ is dumb for some reason)
	numMultiSpreaders not implemented yet, keep at 0. Hopefully will later describe the number of people who interact with all people in a certain residence hall (RA's, dorm heads, etc.)
*/
void create_res_hall(std::vector<int> community_sizes, std::vector<double> community_spread_rates, int numMultiSpreaders, std::string path)
{
	int numcommunity = community_spread_rates.size();
	arma::Mat<float> dir_graph(community_sizes.at(0) + numMultiSpreaders, community_sizes.at(0) + numMultiSpreaders, arma::fill::zeros);
	for (int i = 0; i < numcommunity; i++)
	{
		for (int j = 0; j < community_sizes.at(0)/community_sizes.at(i); j++)
		{
			for (int k = j * community_sizes.at(i); k < (j+1) * community_sizes.at(i) - 1; k++)
			{
				for (int p = k + 1; p < (j + 1) * community_sizes.at(i) - 0; p++)
				{
					dir_graph(k, p) = dir_graph(p, k) = community_spread_rates.at(i);
				}
			}
		}
	}
	//To implement for multi-spreader code later
	/*for (int q = 0; q < numSuperSpreaders; q++)
	{

	}
	*/
	dir_graph.save(path);
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
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	int currentCommunitySize;
	int currentCommunityWeight;
	int pop = mat.n_cols;
	std::string printhead = "Population is: ";
	std::cout << printhead + std::to_string(pop);
	std::vector<int> randarray(pop);
	for (int j = 0; j < pop; j++)
	{
		randarray[j] = j;
	}
	for (int i = 0; i < numcommunities; i++)
	{
		int rand_size = int(communitysize_dist(generator));
		float rand_weight = communityweight_dist(generator);
		currentCommunitySize = rand_size < 2 ? 2 : rand_size > pop ? pop : rand_size;
		currentCommunityWeight = rand_weight < 0 ? 0 : rand_weight > 1 ? 1 : rand_weight;
		std::shuffle(std::begin(randarray), std::end(randarray), generator);

		for (int k = 0; k < currentCommunitySize; k++) 
		{
			for (int m = k + 1; m < currentCommunitySize; m++) 
			{

				mat(randarray[k], randarray[m]) = mat(randarray[m], randarray[k]) = mat(randarray[k], randarray[m]) + currentCommunityWeight > 1 ? 1 : mat(randarray[k], randarray[m]) + currentCommunityWeight;
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
	std::string head = "\nConcatenating ";
	for (fsys::directory_iterator mat{ path }; mat != end; mat++) 
	{
		std::cout << (head + mat->path().string());
		cache.load(mat->path().string(), arma::auto_detect);
		int init_dim = concat_mat.n_cols;
		concat_mat.resize(init_dim + cache.n_cols, init_dim + cache.n_cols);
		for (int i = init_dim; i < init_dim + cache.n_cols; i++)
		{
			for (int j = i; j < init_dim + cache.n_cols; j++) 
			{
				concat_mat(i, j) = concat_mat(j, i) = cache(i - init_dim, j - init_dim);
			}
		}
		
	}
	for (int k = 0; k < concat_mat.n_cols; k++) 
	{
		concat_mat(k, k) = 1;
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
