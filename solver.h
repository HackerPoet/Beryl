#pragma once
#include "util.h"
#include <random>

//using RNG = std::mt19937;
using RNG = std::minstd_rand;
extern RNG eng;

void set_rand_seed(int seed);
void add_random_noise(const VectorXf& x, float sigma, VectorXf& result);
void set_random_mags(VectorXf& x, float sigma);
void initialize_random_noise(VectorXf& x, float sigma, size_t size);
void initialize_random_noise(std::vector<std::pair<float, VectorXf>>& xvec, float sigma);
void initialize_random_noise(std::vector<std::pair<float, VectorXf>>& xvec, const VectorXf& guess, float sigma);

float main_optimizer(float (*objective_function)(const Verts3D&, const Planes&), VectorXf& result, int max_iters,
                     float sigma, float beta, int clusters);
bool study_sample(float (*objective_function)(const Verts3D&, const Planes&), VectorXf& result, int max_iters,
                  int clusters, float sigma, float beta, bool aggressive, float smallintMul, float q_score);
bool explore_shape(const std::string& load_fname, const std::string& save_dir, int max_iters, float smallintMul);
