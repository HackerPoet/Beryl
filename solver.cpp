#include "solver.h"
#include "globals.h"
#include <algorithm>
#include <iostream>
#include <filesystem>

RNG eng;

struct Sample {
    VectorXf x;
    float score;
    int time;

    bool operator<(const Sample& rhs) const {
        return this->score < rhs.score;
    }
};

void set_rand_seed(int seed) {
  eng.seed(seed);
}

void add_random_noise(const VectorXf& x, float sigma, VectorXf& result) {
    std::normal_distribution<float> rand_normal(0.0f, sigma);
    result = x;
    const int num_k = std::uniform_int_distribution<int>(1, 3)(eng);
    for (int k = 0; k < num_k; ++k) {
        if (g_sym.perms.empty()) {
            const int i = std::uniform_int_distribution<int>(0, (int)x.size() / 3 - 1)(eng);
            result[i * 3 + 0] += rand_normal(eng);
            result[i * 3 + 1] += rand_normal(eng);
            result[i * 3 + 2] += rand_normal(eng);
        } else {
            const Face& group = g_sym.perms[std::uniform_int_distribution<size_t>(0, g_sym.perms.size() - 1)(eng)].cycle;
            for (int i : group) {
                result[i * 3 + 0] += rand_normal(eng);
                result[i * 3 + 1] += rand_normal(eng);
                result[i * 3 + 2] += rand_normal(eng);
            }
        }
    }
}

void set_random_mags(VectorXf& x, float sigma) {
  std::normal_distribution<float> rand_normal(0.0f, sigma);
  for (int i = 0; i < x.size(); i += 3) {
    Eigen::Map<Vector3f> sub_x(x.data() + i);
    sub_x.normalize();
    sub_x *= rand_normal(eng);
  }
}

void initialize_random_noise(VectorXf& x, float sigma, size_t size) {
    std::normal_distribution<float> rand_normal(0.0f, sigma);
    x.resize(size);
    for (int i = 0; i < x.size(); ++i) {
        x[i] = rand_normal(eng);
    }
}

void initialize_random_noise(std::vector<std::pair<float, VectorXf>>& xvec, float sigma) {
    size_t size = (g_dual ? g_tris.size() : g_polys.size()) * 3;
    for (auto& x : xvec) {
        initialize_random_noise(x.second, sigma, size);
    }
}

void initialize_random_noise(std::vector<std::pair<float, VectorXf>>& xvec, const VectorXf& guess, float sigma) {
    for (auto& x : xvec) {
        initialize_random_noise(x.second, sigma, guess.size());
        x.second += guess;
    }
}

void initialize_random_mag(std::vector<VectorXf>& yvec, const VectorXf& n, float sigma) {
    std::normal_distribution<float> rand_normal(0.0f, sigma);
    for (VectorXf& y : yvec) {
        y.resize(g_polys.size());
        for (int i = 0; i < y.size(); ++i) {
            y[i] = rand_normal(eng);
        }
    }
}

bool v_pred(const std::pair<float, VectorXf>& left, const std::pair<float, VectorXf>& right) {
    return left.first < right.first;
}

float main_optimizer(float (*objective_function)(const Verts3D&, const Planes&), VectorXf& result, int max_iters,
                   float sigma, float beta, int clusters) {
    static const int extra_tries = 1;
    static Planes planes;
    static Verts3D v3ds;
    std::vector<std::pair<float, VectorXf>> xv(clusters * extra_tries);
    VectorXf new_pt;

    initialize_random_noise(xv, sigma);
    for (auto& x : xv) {
      if (!g_sym.Apply(x.second)) {
        return -1.0f;
      }
      if (g_dual) {
        y_to_v3ds(x.second, v3ds);
        v3ds_to_planes(v3ds, g_polys, planes);
      } else {
        x_to_planes(x.second, planes);
        planes_to_v3ds(g_tris, planes, v3ds);
      }
      x.first = objective_function(v3ds, planes);
    }
    std::sort(xv.begin(), xv.end(), v_pred);
    xv.resize(clusters);

    int iter = 0;
    int last_best_iter = 0;
    float best_cost = 99999;
    const size_t x_size = xv[0].second.size();
    new_pt.resize(x_size);
    while (true) {
        const size_t min_ix = std::min_element(xv.begin(), xv.end(), v_pred) - xv.begin();
        const float min_cost = xv[min_ix].first;
        iter += 1;
        if (min_cost < best_cost) {
            std::cout << min_cost << "    " << iter << std::endl;
            best_cost = min_cost;
            result = xv[min_ix].second;
            last_best_iter = iter;
        }

        if (iter - last_best_iter > max_iters) { break; }
        if (best_cost == 0.0f) { break; }

        for (size_t i = 0; i < clusters; ++i) {
            static const float gamma = 1.25f;
            size_t ix = i;
            if (xv[i].first > best_cost * gamma) {
                ix = std::uniform_int_distribution<size_t>(0, clusters-1)(eng);
            }
            add_random_noise(xv[ix].second, 1.0f, new_pt);
            new_pt *= beta;
            if (!g_sym.Apply(new_pt)) {
              return -1.0f;
            }
            if (g_dual) {
              y_to_v3ds(new_pt, v3ds);
              v3ds_to_planes(v3ds, g_polys, planes);
            } else {
              x_to_planes(new_pt, planes);
              planes_to_v3ds(g_tris, planes, v3ds);
              if (!is_finite(v3ds)) { continue; }
            }
            const float new_cost = objective_function(v3ds, planes);
            float cost_mult = 1.25f; //1.0f + 0.05f * std::uniform_int_distribution<int>(0, 5)(eng);
            if (i == min_ix) { cost_mult = 1.0f; }
            if ((new_cost <= xv[i].first) || (new_cost <= best_cost * cost_mult)) {
              xv[i].first = new_cost;
              xv[i].second = new_pt;
            }
        }
    }
    return best_cost;
}

void truncate_normalize(VectorXf& result, float smallintMul) {
  if (smallintMul > 0.0f) {
    result *= (smallintMul * result.size() / result.norm());
    trunc_x(result);
  } else {
    result *= (result.size() / result.norm());
  }
}

bool study_sample(float (*objective_function)(const Verts3D&, const Planes&), VectorXf& result, int max_iters, int clusters, float sigma, float beta, bool aggressive, float smallintMul, float q_score) {
    static Planes planes;
    static Verts3D v3ds;  
  
    const size_t x_size = result.size();
    if (!g_sym.Apply(result)) { return false; }
    truncate_normalize(result, smallintMul);
    if (g_dual) {
        y_to_v3ds(result, v3ds);
        v3ds_to_planes(v3ds, g_polys, planes);
        if (!is_finite(planes)) {
          std::cout << "ERROR: NaN values in polyhedron." << std::endl;
          return false;
        }
    } else {
        x_to_planes(result, planes);
        planes_to_v3ds(g_tris, planes, v3ds);
        if (!is_finite(v3ds)) {
          std::cout << "ERROR: NaN values in polyhedron." << std::endl;
          return false;
        }
    }

    float min_score = objective_function(v3ds, planes) + q_penalty(v3ds, planes, count_crossings(v3ds, planes));
    std::vector<VectorXf> xv(clusters, result);
    std::vector<float> scores(clusters, min_score);
    VectorXf new_x(x_size);

    for (int iter = 0; iter < max_iters; ++iter) {
        int num_updated = 0;
        for (int i = 0; i < clusters; ++i) {
            add_random_noise(xv[i], sigma, new_x);
            if (!g_sym.Apply(new_x)) { return false; }
            truncate_normalize(new_x, smallintMul);
            if (g_dual) {
                y_to_v3ds(new_x, v3ds);
                v3ds_to_planes(v3ds, g_polys, planes);
                if (!is_finite(planes)) { continue; }
            } else {
                x_to_planes(new_x, planes);
                planes_to_v3ds(g_tris, planes, v3ds);
                if (!is_finite(v3ds)) { continue; }
            }
            float new_cost = objective_function(v3ds, planes);
            if (new_cost < scores[i] || (!aggressive && int(new_cost) <= int(scores[i]))) {
                //Only compute costly q_penalty when score is low enough to matter
                if (int(new_cost) <= int(q_score)) {
                    new_cost += q_penalty(v3ds, planes, count_crossings(v3ds, planes));
                }
                if (new_cost < min_score) {
                    std::cout << "==== New Best! ==== (" << new_cost << ")" << std::endl;
                    result = new_x;
                    if (int(new_cost) < int(min_score) || new_cost <= 1.0f) {
                        std::fill(xv.begin(), xv.end(), new_x);
                        std::fill(scores.begin(), scores.end(), new_cost);
                    }
                    min_score = new_cost;
                }
                if (!aggressive || new_cost == min_score) {
                    xv[i] = new_x;
                    scores[i] = new_cost;
                    num_updated += 1;
                }
            }
        }
        if (iter % 10 == 0) {
            std::cout << "[" << iter << "] Updated: " << num_updated << "/" << clusters << "      " << sigma << std::endl;
        }
        if (num_updated < clusters / 10) {
            sigma *= beta;
        } else if (num_updated > clusters / 5) {
            sigma *= 1.005f;
        }
        if (sigma < 5e-7f) { break; }
    }
    return true;
}

bool explore_shape(const std::string& load_fname, const std::string& save_dir, int max_iters, float smallintMul) {
  //Import an example obj file
  int f_iter = 0;
  Verts3D obj_verts;
  Planes obj_planes;
  VectorXf obj_x;
  if (load_fname.empty()) {
    dual_graph(g_tris, g_polys, g_edges);
    fix_face_ordering(g_polys, g_edges);
    if (g_dual) {
      std::swap(g_tris, g_polys);
      make_edges(g_polys, g_edges);
      fix_face_ordering(g_polys, g_edges);
    }
    initialize_random_noise(obj_x, 1.0f, g_polys.size() * 3);
    if (g_dual) {
      y_to_v3ds(obj_x, obj_verts);
    } else {
      x_to_planes(obj_x, obj_planes);
      planes_to_v3ds(g_tris, obj_planes, obj_verts);
    }
  } else {
    std::vector<std::string> name_split = split(std::filesystem::path(load_fname).stem().string(), '_');
    f_iter = std::atoi(name_split[std::min(size_t(3), name_split.size() - 1)].c_str());
    if (!import_obj(load_fname.c_str(), obj_verts, g_polys)) {
      std::cout << "ERROR: Failed to load obj file: " << load_fname << std::endl;
      return false;
    }
  }
#if 1
  Edges dual_edges;
  make_edges(g_polys, g_edges);
  dual_graph(g_polys, g_tris, dual_edges);
  v3ds_to_planes(obj_verts, g_polys, obj_planes);
  if (g_dual) {
    v3ds_to_y(obj_verts, obj_x);
  } else {
    planes_to_x(obj_planes, obj_x);
  }
#else
  dual_graph(g_polys, g_tris, g_edges);
  make_edges(g_polys, g_edges);
  std::swap(g_polys, g_tris);
  v3ds_to_y(obj_verts, obj_x);
  x_to_planes(obj_x, obj_planes);
  planes_to_v3ds(g_tris, obj_planes, obj_verts);
#endif

  //Print characteristics
  std::cout << "Petrie Length: " << petrie_length(g_dual ? g_polys : g_tris) << std::endl;
  std::cout << "===================" << std::endl;
  if (!load_fname.empty()) {
    std::cout << "Loaded:        " << load_fname << std::endl;
  }
  save_sample("refined", obj_planes, obj_verts, f_iter, false);
  std::cout << "===================" << std::endl;
  auto objective = (g_dual ? objective_int : objective_sum);
  float best_score = objective(obj_verts, obj_planes) + q_penalty(obj_verts, obj_planes, count_crossings(obj_verts, obj_planes));
  float noise_level = 8.0f;

  while (true) {
    VectorXf new_x = obj_x;
    new_x *= new_x.size() / new_x.norm();
    if (best_score >= 1.0f) {
      std::cout << "Noise level: " << noise_level << std::endl;
      std::normal_distribution<float> rand_noise(0.0f, 1.0f);
      for (int i = 0; i < new_x.size() / 3; ++i) {
        if (rand_noise(eng) > 1.25f) {
          new_x[i * 3 + 0] += noise_level * rand_noise(eng);
          new_x[i * 3 + 1] += noise_level * rand_noise(eng);
          new_x[i * 3 + 2] += noise_level * rand_noise(eng);
        } else {
          new_x[i * 3 + 0] += noise_level * 0.05f * rand_noise(eng);
          new_x[i * 3 + 1] += noise_level * 0.05f * rand_noise(eng);
          new_x[i * 3 + 2] += noise_level * 0.05f * rand_noise(eng);
        }
      }
    }

    //Run optimizer
    if (!study_sample(objective, new_x, max_iters, 100, 10.0f, 0.997f, best_score < 1.0f, smallintMul, best_score)) {
      return false;
    }
    if (g_dual) {
      y_to_v3ds(new_x, obj_verts);
      v3ds_to_planes(obj_verts, g_polys, obj_planes);
    } else {
      x_to_planes(new_x, obj_planes);
      planes_to_v3ds(g_tris, obj_planes, obj_verts);
    }
    const float new_score = objective(obj_verts, obj_planes) + q_penalty(obj_verts, obj_planes, count_crossings(obj_verts, obj_planes));
    if (int(new_score) <= int(best_score)) {
      noise_level *= 1.5f;
    } else {
      noise_level *= 0.85f;
    }
    if (new_score < best_score) {
      std::cout << "########################################" << std::endl;
      std::cout << "  NEW BEST SCORE: " << best_score << " --> " << new_score << std::endl;
      std::cout << "########################################" << std::endl;
      obj_x = new_x;
      best_score = new_score;
    } else {
      std::cout << "Best score: " << best_score << std::endl;
    }
    save_sample((save_dir + "/refined").c_str(), obj_planes, obj_verts, f_iter, best_score < 1.0f || new_score == best_score);
  }
  return true;
}
