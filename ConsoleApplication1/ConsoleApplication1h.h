#pragma once
#include <armadillo>

class cascadingmodel
{
public:
	void onetime_load();
	arma::Col<int> step(arma::Col<int> currentInfected);
	void create_res_hall(int community_sizes[], double community_spread_rates[], int numcommunity, int numMultiSpreaders, std::string path);
	arma::Mat<float> create_communities(arma::Mat<float> mat, int numcommunities, int avg_communitysize, float communitysize_variation, float communityweight, float communityweight_variation);
	arma::Mat<float> concat_diag_mat(std::experimental::filesystem::path path);
};