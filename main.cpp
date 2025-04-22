#include "cmdLine.h"
#include "util.h"
#include "solver.h"
#include "globals.h"
#include <chrono>
#include <iostream>
#include <fstream>
#include <filesystem>

//Global command line arguments and default values
bool        h_printHelp = false;
bool        h_printEmbed = false;
bool        h_printEmbedV = false;
bool        h_printEmbedS = false;
std::string h_map = "";
bool        h_zeroIndex = false;
std::string h_outputDir = "";
bool        h_dual = false;
std::string h_perms = "";
int         h_iters = 1600;
int         h_thresh = 0;
std::string h_load = "";
float       h_smallint = 0.0f;
std::string h_cutout = "";
std::string h_wireframe = "";
bool        h_printsym = false;

//Start of program
int main(int argc, const char* argv[]) {
  CmdLine cmd("Beryl");
  cmd.addCategory("============= General =============");
  cmd.addArgument({"-h","--help"}, "Print this help message", &h_printHelp);
  cmd.addArgument({"-m","--map"}, "Name of map to load triangles and automorphisms", &h_map);
  cmd.addArgument({"-e","--embeddings"}, "Print possible embeddable geometric symmetries", &h_printEmbed);
  cmd.addArgument({"-ev","--verbose-embeddings"}, "Print all geometric symmetries including intersections", &h_printEmbedV);
  cmd.addArgument({"-es","--save-embeddings"}, "Save all geometric symmetries into the Map folder", &h_printEmbedS);
  cmd.addArgument({"-z","--zindex"}, "Map files are zero-indexed if set", &h_zeroIndex);

  cmd.addCategory("============= Solver =============");
  cmd.addArgument({"-o","--output"}, "Directory to save OBJ files", &h_outputDir);
  cmd.addArgument({"-d","--dual"}, "Solve for dual embedding instead of triangles", &h_dual);
  cmd.addArgument({"-s","--sym"}, "Geometric symmetry to use and list of permutations", &h_perms);
  cmd.addArgument({"-i","--iters"}, "End solve if no improvement after this many iterations", &h_iters);
  cmd.addArgument({"-t","--thresh"}, "Only save if this many intersections or fewer", &h_thresh);

  cmd.addCategory("============= Refinement =============");
  cmd.addArgument({"-l","--load"}, "OBJ file to load", &h_load);
  cmd.addArgument({"-o","--output"}, "Directory to save OBJ files", &h_outputDir);
  cmd.addArgument({"-i","--iters" }, "Fixed number of iterations per refinement", &h_iters);
  cmd.addArgument({"-si","--smallint"}, "Use small integer coordinates, multiply by this amount", &h_smallint);

  cmd.addCategory("============= Other =============");
  cmd.addArgument({"-l","--load"}, "Input OBJ file to load", &h_load);
  cmd.addArgument({"-c","--cutout"}, "Output filename for OBJ cutout", &h_cutout);
  cmd.addArgument({"-w","--wireframe"}, "Output filename for OBJ wireframe", &h_wireframe);
  cmd.addArgument({"-p","--printsym"}, "Print embedded symmetries of the OBJ file", &h_printsym);

  //Parse the command line
  if (!cmd.parse(argc, argv)) {
    return 1;
  }

  //Setup filenames
  g_dual = !h_dual;
  std::string triFile = "Maps/" + h_map + ".txt";
  std::string amFile = "Maps/" + h_map + "_automorphisms.txt";
  SymmetryList syms;

  //Setup output directory automatically if none specified
  if (h_outputDir.empty()) {
    std::string symId = h_perms;
    std::replace(symId.begin(), symId.end(), ',', '_');
    h_outputDir = (g_dual ? "Tri/" : "Poly/") + h_map + "/" + symId;
  }

  if (h_printHelp || argc <= 1) {
    //Print the help menu
    cmd.printHelp();
    return 0;
  } else if (h_printEmbedS) {
    static const std::string suffix("_automorphisms.txt");
    for (const auto& dir_entry : std::filesystem::directory_iterator("Maps")) {
      //Load each regular map and automorphism file
      const auto& path = dir_entry.path();
      const std::string amName = path.filename().string();
      if (ends_with(amName, suffix)) {
        h_map = amName.substr(0, amName.length() - suffix.length());
        triFile = "Maps/" + h_map + ".txt";
        amFile = "Maps/" + h_map + "_automorphisms.txt";
        const std::string outFile = "Maps/" + h_map + "_symmetries.txt";
        if (!open_face_file(triFile.c_str(), g_tris, h_zeroIndex)) {
          std::cout << "ERROR: Failed to load triangle list: " << triFile << std::endl;
          return 1;
        }

        //Print the basic statistics
        dual_graph(g_tris, g_polys, g_edges);
        const int petrieLength = petrie_length(g_tris);
        std::cout << "V,E,F: " << g_polys.size() << "," << g_edges.size() << "," << g_tris.size() << std::endl;
        std::cout << "Schlafi: {" << g_tris[0].size() << "," << g_polys[0].size() << "}_" << petrieLength << std::endl;

        if (!syms.Load(g_tris, amFile.c_str())) {
          std::cout << "ERROR: Failed to load automorphism file: " << amFile << std::endl;
          return 1;
        }
        std::cout << "Saving all unique symmetries: " << h_map << std::endl;
        std::ofstream fout(outFile);
        syms.Analyze(g_tris, true, fout);
      }
    }
    return 0;
  } else if (h_printEmbedV || h_printEmbed) {
    //Make sure a triangle map is specified
    if (h_map.empty()) {
      std::cout << "ERROR: Must specify a triangle list (map)." << std::endl;
      return 1;
    }

    //Open the triangle file
    if (!open_face_file(triFile.c_str(), g_tris, h_zeroIndex)) {
      std::cout << "ERROR: Failed to load triangle list: " << triFile << std::endl;
      return 1;
    }

    //Print the basic statistics
    dual_graph(g_tris, g_polys, g_edges);
    const int petrieLength = petrie_length(g_tris);
    std::cout << "V,E,F: " << g_polys.size() << "," << g_edges.size() << "," << g_tris.size() << std::endl;
    std::cout << "Schlafi: {" << g_tris[0].size() << "," << g_polys[0].size() << "}_" << petrieLength << std::endl;

    //Open the automorphism file
    if (!syms.Load(g_tris, amFile.c_str())) {
      std::cout << "ERROR: Failed to load automorphism file: " << amFile << std::endl;
      return 1;
    }

    //Run the analysis
    std::cout << "Printing all unique symmetries..." << std::endl;
    syms.Analyze(g_tris, !h_printEmbedV);
    return 0;
  } else if (!h_cutout.empty()) {
    //Save a cutout of the shape
    Verts3D vert3ds;
    Planes planes;
    if (!import_obj(h_load.c_str(), vert3ds, g_polys)) {
      std::cout << "ERROR: Failed to load obj file: " << h_load << std::endl;
      return 1;
    }
    v3ds_to_planes(vert3ds, g_polys, planes);
    if (!export_cutout(h_cutout.c_str(), vert3ds, planes, 10.0f)) {
      std::cout << "ERROR: Failed to save cutout obj file: " << h_cutout << std::endl;
      return 1;
    }
    return 0;
  } else if (!h_wireframe.empty()) {
    //Save a cutout of the shape
    Verts3D vert3ds;
    if (!import_obj(h_load.c_str(), vert3ds, g_polys)) {
      std::cout << "ERROR: Failed to load obj file: " << h_load << std::endl;
      return 1;
    }
    if (!export_wireframe_obj(h_wireframe.c_str(), vert3ds, g_polys, 1.0f)) {
      std::cout << "ERROR: Failed to save wireframe obj file: " << h_wireframe << std::endl;
      return 1;
    }
    return 0;
  } else if (h_printsym) {
    //Print symmetries of the shape
    Verts3D vert3ds;
    if (!import_obj(h_load.c_str(), vert3ds, g_tris)) {
      std::cout << "ERROR: Failed to load obj file: " << h_load << std::endl;
      return 1;
    }
    dual_graph(g_tris, g_polys, g_edges);
    if (!g_dual) {
      std::swap(g_polys, g_tris);
    }
    if (h_map.empty()) {
      std::cout << "ERROR: Must specify automorphism file." << std::endl;
      return 1;
    }
    if (!syms.Load(g_tris, amFile.c_str())) {
      std::cout << "ERROR: Failed to load automorphism file: " << amFile << std::endl;
      return 1;
    }
    std::cout << "Searching..." << std::endl;
    for (size_t i = 0; i < syms.syms.size(); ++i) {
      const std::string result = syms.syms[i].Test(vert3ds);
      if (!result.empty()) {
        std::cout << result << " " << i << std::endl;
      }
    }
    return 0;
  }

  //Setup seed for the solver
  const auto time_now = std::chrono::high_resolution_clock::now().time_since_epoch();
  const unsigned int seed = (unsigned int)std::chrono::duration_cast<std::chrono::milliseconds>(time_now).count();
  set_rand_seed(seed);

  //Load the triangle list
  if (h_map.empty()) {
    std::cout << "ERROR: Must specify a triangle list (map)." << std::endl;
    return 1;
  } else if (!open_face_file(triFile.c_str(), g_tris, h_zeroIndex)) {
    std::cout << "ERROR: Failed to load triangle list: " << triFile << std::endl;
    return 1;
  }

  //Load symmetry
  if (h_perms.length() == 0) {
    std::cout << "WARNING: No permutation groups specified. Solver will not use any symmetry." << std::endl;
  } else {
    const SymType type = Symmetry::GetType(h_perms[0]);
    if (type == SymType::None) {
      std::cout << "ERROR: Unknown symmetry type: " << h_perms[0] << std::endl;
      return 1;
    }
    if (h_perms.length() == 1) {
      std::cout << "ERROR: No symmetry index list specified." << std::endl;
      return 1;
    }
    const int prodSize = Symmetry::GetProdSize(type);
    std::vector<std::string> permStrs = split(h_perms.substr(1, h_perms.size() - 1), ',');
    if (prodSize != (int)permStrs.size()) {
      std::cout << "ERROR: Symmetry product count doesn't match type '" << h_perms[0] << "'." << std::endl;
      std::cout << "       Expected " << prodSize << " group(s) but given " << permStrs.size() << "." << std::endl;
      return 1;
    }
    std::vector<int> permIxs;
    for (const std::string& s : permStrs) {
      try {
        permIxs.push_back(std::stoi(s));
      } catch (...) {
        std::cout << "ERROR: Could not parse comma-separated permutation indices: " << h_perms << std::endl;
        return 1;
      }
    }
    if (!Symmetry::Load(g_sym, g_sym_tri, g_tris, amFile.c_str(), permIxs, h_zeroIndex)) {
      std::cout << "ERROR: Failed to load automorphism file: " << amFile << std::endl;
      return 1;
    }
    g_sym.SetType(type);
  }

  //Check if this is a refinement solve
  if (!h_load.empty()) {
    std::cout << "Max Iters: " << h_iters << std::endl;
    if (!explore_shape(h_load, h_outputDir, h_iters, h_smallint)) {
      std::cout << "ERROR: Could not refine polyhedron." << std::endl;
    }
    return 1;
  }

  //Find the dual graph to get the polygon and edge linkage
  dual_graph(g_tris, g_polys, g_edges);
  fix_face_ordering(g_polys, g_edges);
  if (g_dual) {
    std::swap(g_tris, g_polys);
    make_edges(g_polys, g_edges);
    fix_face_ordering(g_polys, g_edges);
  }

  //Print information for debugging purposes
  const int petrieLength = petrie_length(g_dual ? g_polys : g_tris);
  std::cout << "V,E,F: " << g_tris.size() << "," << g_edges.size() << "," << g_polys.size() << std::endl;
  std::cout << "Schlafi: {" << g_polys[0].size() << "," << g_tris[0].size() << "}_" << petrieLength << std::endl;
  std::cout << "Max Seq Iters: " << h_iters << std::endl;

  //Start the main solver
  int iter = 0;
  int best_score = 99999999;
  while (true) {
    //Create directory for results
    if (!std::filesystem::exists(h_outputDir)) {
      std::filesystem::create_directories(h_outputDir);
    }

    //Run the optimizer
    VectorXf result;
    float score = main_optimizer((g_dual ? objective_int_qlim : objective_sum), result, h_iters, 1.0f, 0.996f, 16);
    if (score < 0.0f) {
      std::cout << "ERROR: Solver failed to run." << std::endl;
      return 1;
    }

    //Get the actual values of crossings and intersection independent of score
    Planes planes;
    Verts3D v3ds;
    if (g_dual) {
      y_to_v3ds(result, v3ds);
      v3ds_to_planes(v3ds, g_polys, planes);
    } else {
      x_to_planes(result, planes);
      planes_to_v3ds(g_tris, planes, v3ds);
    }
    const int crossings = count_crossings(v3ds, planes);
    const int intersections = count_intersections(v3ds, planes);
    const int sum = crossings + intersections;

    //Check we should save it
    std::cout << "Score    : " << score << std::endl;
    if (is_finite(v3ds)) {
      const bool can_save = (sum <= h_thresh || (h_thresh <= 0 && sum <= best_score));
      const std::string save_str = h_outputDir + "/shape";
      save_sample(save_str.c_str(), planes, v3ds, iter, can_save);
    }
    iter += 1;

    //Quit once a solution has been found
    if (crossings + intersections == 0) { break; }
    if (sum < best_score) { best_score = sum; }
  }
  return 0;
}
