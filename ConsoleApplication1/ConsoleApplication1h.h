#pragma once
#include <armadillo>
#include <map>
#include <vector>

class cascadingmodel
{
public:
	arma::Col<int> step(arma::Col<int> currentInfected, arma::Col<int> days_to_removal, int max_days_to_removal);
	void onetime_load();
	void create_res_hall(std::vector<int> community_sizes, std::vector<double> community_spread_rates, int numMultiSpreaders, std::string path);
	arma::Mat<float> create_communities(arma::Mat<float> mat, int numcommunities, int avg_communitysize, float communitysize_variation, float communityweight, float communityweight_variation);
	arma::Mat<float> concat_diag_mat(std::experimental::filesystem::path path);
	arma::Mat<float> remove_nodes(arma::Mat<float> initial_mat, arma::Col<int> days_to_removal);
	void create_dir_graph(std::map<std::string, float> community_info, std::map<std::string, std::vector<int>> apts, std::map<int, std::vector<double>> apts_rates, std::map<std::string, std::vector<int>> suite, std::map<int, std::vector<double>> suite_rates, std::map<std::string, std::vector<int>>trad, std::map<int, std::vector<double>> trad_rates, std::map<std::string, std::vector<int>> greek, std::map<int, std::vector<double>> greek_rates, std::string path);

};