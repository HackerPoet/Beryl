#include "util.h"
#include "globals.h"
#include <set>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cfloat>

template<class VTYPE>
inline float pt_to_line_dist_sq(const VTYPE& pt, const VTYPE& a, const VTYPE& b) {
  const VTYPE p = pt - a;
  const VTYPE d = b - a;
  const float t = p.dot(d) / d.squaredNorm();
  const VTYPE v = p - d * std::min(std::max(t, 0.0f), 1.0f);
  return v.squaredNorm();
}

void ltrim(std::string& s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
    return !std::isspace(ch);
  }));
}

void rtrim(std::string& s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
    return !std::isspace(ch);
  }).base(), s.end());
}

bool ends_with(const std::string& value, const std::string& ending) {
  if (ending.size() > value.size()) { return false; }
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

int triangle_id(int s1, int s2, int s3) {
  const int m1 = std::max(s1, std::max(s2, s3));
  const int m3 = std::min(s1, std::min(s2, s3));
  const int m2 = s1 + s2 + s3 - m1 - m3;
  return m1 + 1000 * (m2 + 1000 * m3);
}
int triangle_id(const Face& tri) {
  return triangle_id(tri[0], tri[1], tri[2]);
}

float softmax(float a, float b) {
    return 1.0f - std::exp((std::log((1.0f - a) + 1e-12f) +
                            std::log((1.0f - b) + 1e-12f)) / 2.0f);
}
float softmax(float a, float b, float c) {
    return 1.0f - std::exp((std::log((1.0f - a) + 1e-12f) +
                            std::log((1.0f - b) + 1e-12f) +
                            std::log((1.0f - c) + 1e-12f)) / 3.0f);
}
float softmax(float a, float b, float c, float d) {
    return 1.0f - std::exp((std::log((1.0f - a) + 1e-12f) +
                            std::log((1.0f - b) + 1e-12f) +
                            std::log((1.0f - c) + 1e-12f) +
                            std::log((1.0f - d) + 1e-12f)) / 3.0f);
}

std::vector<std::string> split(const std::string& s, char delim) {
  std::vector<std::string> result;
  std::stringstream ss(s);
  std::string item;
  while (getline(ss, item, delim)) {
    result.push_back(item);
  }
  return result;
}

bool is_finite(const Verts3D& verts) {
    for (const Vector3f& v : verts) {
        if (!v.allFinite()) {
            return false;
        }
    }
    return true;
}
bool is_finite(const Planes& planes) {
    for (const Plane& plane : planes) {
        if (!plane.n.allFinite()) {
            return false;
        }
    }
    return true;
}

void print_faces(const Faces& faces) {
  for (const Face& face : faces) {
    std::cout << "[";
    for (size_t ix = 0; ix < face.size(); ++ix) {
      if (ix > 0) { std::cout << " "; }
      std::cout << face[ix];
    }
    std::cout << "] ";
  }
  std::cout << std::endl;
}

bool open_face_file(const char* fname, Faces& tris, bool zero_indexed) {
    std::cout << "Loading OBJ file: " << fname << "..." << std::endl;
    std::ifstream fin(fname);
    if (!fin.is_open()) { return false; }
    tris.clear();
    std::string line;
    const int zi_sub = (zero_indexed ? 0 : 1);
    while (std::getline(fin, line)) {
        Face face;
        ltrim(line);
        rtrim(line);
        std::vector<std::string> strs = split(line, ' ');
        for (const std::string& str : strs) {
            if (str.size() == 0) { continue; }
            face.push_back(std::stoi(str) - zi_sub);
        }
        tris.push_back(face);
    }
    std::cout << "Loaded " << tris.size() << " triangles." << std::endl;
    return true;
}

bool open_cyclic_file(const char* fname, Faces& tris, bool zero_indexed) {
    std::ifstream fin(fname);
    if (!fin.is_open()) { return false; }
    tris.clear();
    std::string line;
    const int zi_sub = (zero_indexed ? 0 : 1);
    while (std::getline(fin, line)) {
        ltrim(line);
        rtrim(line);
        std::vector<std::string> strs = split(line, ' ');
        const int f1 = std::stoi(strs[0]) - zi_sub;
        for (size_t ix = 1; ix < strs.size(); ++ix) {
            const size_t prev_ix = (ix == 1 ? strs.size() - 1 : ix - 1);
            const int f2 = std::stoi(strs[ix]) - zi_sub;
            const int f3 = std::stoi(strs[prev_ix]) - zi_sub;
            if (f1 < f2 && f1 < f3) {
                tris.push_back(Face{ f1, f2, f3 });
            }
        }
    }
    std::cout << "Loaded " << tris.size() << " triangles." << std::endl;
    return true;
}

inline int alphanumeric_to_int(char c) {
    if (c >= '0' && c <= '9') { return c - '0'; }
    if (c >= 'a' && c <= 'z') { return c - 'a' + 10; }
    if (c >= 'A' && c <= 'Z') { return c - 'a' + 10; }
    return -1;
}

bool open_topology(const char* fname, Faces& tris, int ix) {
    std::ifstream fin(fname);
    if (!fin.is_open()) { return false; }
    tris.clear();
    std::string line;
    for (int i = 0; i <= ix; ++i) {
        std::getline(fin, line);
    }
    std::vector<std::string> strs = split(line, ' ');
    for (size_t i = 1; i < strs.size(); ++i) {
        const std::string& face_txt = strs[i];
        const int f1 = alphanumeric_to_int(face_txt[0]);
        const int f2 = alphanumeric_to_int(face_txt[1]);
        const int f3 = alphanumeric_to_int(face_txt[2]);
        tris.push_back({f1, f2, f3});
    }
    return true;
}

bool verify_topology(const Faces& faces) {
    std::map<int, int> vertex_count;
    std::map<Edge, int> edge_count;
    for (const Face& face : faces) {
        for (size_t i = 0; i < face.size(); ++i) {
            const int v1 = face[i];
            const int v2 = face[(i + 1) % face.size()];
            vertex_count[v1] += 1;
            const int min_v = std::min(v1, v2);
            const int max_v = std::max(v1, v2);
            edge_count[Edge(min_v, max_v)] += 1;
        }
    }
    if (edge_count.size() != faces[0].size() * faces.size() / 2) {
        std::cout << " Edges found: " << edge_count.size() << std::endl;
        std::cout << "But expected: " << (faces[0].size() * faces.size() / 2) << std::endl;
        return false;
    }
    bool success = true;
    for (auto const& x : edge_count) {
        if (x.second != 2) {
            std::cout << " Expected edge to appear twice but edge [" <<
                         x.first.first << "," << x.first.second << "] occured " <<
                         x.second << " times." << std::endl;
            success = false;
        }
    }
    const size_t expected_vcount = (faces[0].size() == 3 ? vertex_count.size() - 1 : 3);
    for (auto const& x : vertex_count) {
        if (x.second != expected_vcount) {
            std::cout << " Vertex[" << x.first << "] expected to appear " << expected_vcount <<
                         " times but found " << x.second << std::endl;
            success = false;
        }
    }
    return success;
}

void save_sample(const char* name, const Planes& planes, const Verts3D& verts, int iter, bool can_save) {
    const int num_crossings = count_crossings(verts, planes);
    const int num_intersections = count_intersections(verts, planes);
    std::cout << "Iter     : " << iter << std::endl;
    std::cout << "Crossing : " << num_crossings << std::endl;
    std::cout << "Intersect: " << num_intersections << std::endl;
    std::cout << "A-Factor : " << angle_penalty(verts) << std::endl;
    std::cout << "D-Factor : " << dist_penalty(verts) << std::endl;
    std::cout << "L-Factor : " << length_penalty(verts) << std::endl;
    std::cout << "P-Factor : " << plane_penalty(planes) << std::endl;
    std::cout << "T-Factor : " << triangle_penalty(verts) << std::endl;
    if (can_save) {
        std::cout << "============ Exporting... ============" << std::endl;
        std::stringstream ss;
        ss << name << "_c" << int(num_crossings) << +"_i" << num_intersections << "_" << iter << ".obj";
        export_obj(ss.str().c_str(), verts, g_polys);
    }
}

void save_dot_graph(const char* fname, const Edges& edges) {
    std::ofstream fout(fname);
    fout << "graph nodes {\n";
    for (const Edge& edge : edges) {
        fout << "  N" << edge.first << " -- N" << edge.second << ";\n";
    }
    fout << "}\n";
}

void dual_graph(const Faces& faces, Faces& dual_faces, Edges& dual_edges) {
    FaceMap verts;
    EdgeMap edges;
    for (size_t fIx = 0; fIx < faces.size(); ++fIx) {
        const Face& face = faces[fIx];
        for (size_t vIx = 0; vIx < face.size(); ++vIx) {
            const int v1 = face[vIx];
            const int v2 = face[(vIx + 1) % face.size()];
            verts[v1].push_back((int)fIx);
            const Edge edge(std::min(v1, v2), std::max(v1, v2));
            edges[edge].push_back((int)fIx);
        }
    }
    dual_faces.clear();
    for (FaceMap::const_iterator it = verts.begin(); it != verts.end(); ++it) {
      dual_faces.push_back(it->second);
    }
    dual_edges.clear();
    for (EdgeMap::const_iterator it = edges.begin(); it != edges.end(); ++it) {
      dual_edges.push_back(Edge(it->second[0], it->second[1]));
    }
}

void dual_verts(const Faces& faces, const Verts3D& verts, Verts3D& dual_verts) {
    dual_verts.clear();
    for (const Face& face : faces) {
        Vector3f v = Vector3f::Zero();
        for (int i : face) {
            v += verts[i];
        }
        dual_verts.push_back(v / (float)face.size());
    }
}

void make_edges(const Faces& faces, Edges& edges) {
    edges.clear();
    for (size_t fIx = 0; fIx < faces.size(); ++fIx) {
        const Face& face = faces[fIx];
        for (size_t vIx = 0; vIx < face.size(); ++vIx) {
              const int v1 = face[vIx];
              const int v2 = face[(vIx + 1) % face.size()];
              const Edge edge(std::min(v1, v2), std::max(v1, v2));
              if (std::find(edges.begin(), edges.end(), edge) == edges.end()) {
                  edges.push_back(edge);
              }
        }
    }
}

bool test_face_ordering(const Faces& polys) {
    std::set<Edge> edges;
    for (const Face& poly : polys) {
        int prevIx = poly[poly.size() - 1];
        for (int ix : poly) {
            const Edge edge(ix, prevIx);
            if (edges.count(edge) > 0) {
                std::cout << edge.first << "," << edge.second << std::endl;
                return false;
            }
            edges.insert(edge);
            prevIx = ix;
        }
    }
    return true;
}

void fix_face_ordering(Faces& polys, const Edges& edges) {
    Faces fixed_polys;
    for (const Face& poly : polys) {
        Face fixed_poly;
        fixed_poly.push_back(poly[0]);
        std::unordered_set<int> poly_set(poly.begin(), poly.end());
        while (fixed_poly.size() < poly.size()) {
            for (const Edge& e : edges) {
                const int e1 = e.first;
                const int e2 = e.second;
                const int lastVert = fixed_poly[fixed_poly.size() - 1];
                if (e1 == lastVert && poly_set.find(e2) != poly_set.end() &&
                    std::find(fixed_poly.begin(), fixed_poly.end(), e2) == fixed_poly.end()) {
                    fixed_poly.push_back(e2);
                    break;
                }
                if (e2 == lastVert && poly_set.find(e1) != poly_set.end() &&
                    std::find(fixed_poly.begin(), fixed_poly.end(), e1) == fixed_poly.end()) {
                    fixed_poly.push_back(e1);
                    break;
                }
            }
        }
        fixed_polys.push_back(fixed_poly);
    }
    polys = fixed_polys;
}

bool import_obj(const char* fname, Verts3D& verts, Faces& polys) {
    verts.clear();
    polys.clear();
    std::ifstream fin(fname);
    if (!fin.is_open()) { return false; }
    std::string line;
    while (std::getline(fin, line)) {
        ltrim(line);
        rtrim(line);
        std::vector<std::string> strs = split(line, ' ');
        if (strs[0] == "v") {
            const Vector3f vert((float)std::stod(strs[1]), (float)std::stod(strs[2]), (float)std::stod(strs[3]));
            verts.push_back(vert);
        } else if (strs[0] == "f") {
            Face face;
            for (size_t i = 1; i < strs.size(); ++i) {
                const std::string faceStr = strs[i];
                size_t slash_pos = faceStr.find('/', 0);
                if (slash_pos == std::string::npos) {
                    face.push_back(std::stoi(faceStr) - 1);
                } else {
                    face.push_back(std::stoi(faceStr.substr(0, slash_pos)) - 1);
                }
            }
            polys.push_back(face);
        }
    }
    return true;
}

bool export_obj(const char* fname, const Verts3D& verts, const Faces& polys) {
    std::ofstream fout(fname);
    if (!fout.is_open()) { return false; }
    for (Vector3f vert : verts) {
        fout << "v " << std::setprecision(12) << (double)vert.x() << " " << (double)vert.y() << " " << (double)vert.z() << "\n";
    }
    for (Face poly : polys) {
        fout << "f";
        for (int p : poly) {
          fout << " " << (p+1);
        }
        fout << "\n";
    }
    return true;
}

bool export_colored_obj(const char* fname, const Verts3D& verts, const Faces& polys) {
    std::ofstream fout(fname);
    if (!fout.is_open()) { return false; }
    std::map<int, int> tri_map;
    int num_colors = 0;
    const Symmetry& sym = (g_dual ? g_sym_tri : g_sym);
    for (const auto& group : sym.perms) {
        for (int i : group.cycle) {
            tri_map[i] = num_colors;
        }
        num_colors += 1;
    }
    std::cout << "Unique Colors: " << num_colors << std::endl;
    for (size_t i = 0; i < polys.size(); ++i) {
        for (int vix : polys[i]) {
            const Vector3f& vert = verts[vix];
            fout << "v " << std::setprecision(12) << (double)vert.x() << " " << (double)vert.y() << " " << (double)vert.z() << "\n";
        }
    }
    for (size_t i = 0; i < polys.size(); ++i) {
        const float uv = (float)((tri_map[(int)i] * 23) % num_colors);
        for (int j : polys[i]) {
            fout << "vt " << ((uv + 0.5f) / (float)num_colors) << " " << 0.5f << "\n";
        }
    }
    int vert_ix = 1;
    for (size_t i = 0; i < polys.size(); ++i) {
        fout << "f";
        for (int j = 0; j < polys[i].size(); ++j) {
             fout << " " << vert_ix << "/" << vert_ix;
             vert_ix += 1;
        }
        fout << "\n";
    }
    return true;
}

void make_perp_line(const Vector3f& a, const Vector3f& b, const Vector3f& c, Vector3f& pa, Vector3f& pb, float extrude) {
    const Vector3f ba = b - a;
    const Vector3f ca = c - a;
    const Vector3f perp = ba.cross(ca).cross(ba).normalized();
    pa = a + perp * extrude;
    pb = b + perp * extrude;
}

int get_symm_ix(int ix) {
    if (g_sym.perms.size() == 0) { return ix; }
    for (const auto& group : g_sym.perms) {
        if (group.indexOf(ix) >= 0) {
            return group.rep();
        }
    }
    std::cout << "ERROR!" << std::endl;
    return -1;
}

bool export_wireframe_obj(const char* fname, const Verts3D& verts, const Faces& polys, float extrude) {
  std::cout << "Exporting wireframe: " << fname << "..." << std::endl;
  Faces wire_quads;
  Verts3D wire_verts;
  std::map<int, int> face_check;
  std::vector<int> wire_uvs;
  for (const Face& poly : polys) {
    if (poly.size() != 3) {
      std::cout << "ERROR: Wireframe export only supported on triangular meshes." << std::endl;
      return false;
    }
    const int s1 = get_symm_ix(poly[0]);
    const int s2 = get_symm_ix(poly[1]);
    const int s3 = get_symm_ix(poly[2]);
    const int m = triangle_id(s1, s2, s3);
    const auto& existing_face = face_check.find(m);
    if (existing_face == face_check.end()) {
      face_check[m] = (int)face_check.size();
    }
    const int uv = face_check[m];
    const Vector3f& a = verts[poly[0]];
    const Vector3f& b = verts[poly[1]];
    const Vector3f& c = verts[poly[2]];
    const float dsq1 = pt_to_line_dist_sq(a, b, c);
    const float dsq2 = pt_to_line_dist_sq(b, c, a);
    const float dsq3 = pt_to_line_dist_sq(c, a, b);
    float dist = std::sqrt(std::min(std::min(dsq1, dsq2), dsq3));
    dist = std::min(dist * 0.4f, extrude);
    Vector3f pab_a, pab_b;
    Vector3f pbc_b, pbc_c;
    Vector3f pca_c, pca_a;
    make_perp_line(a, b, c, pab_a, pab_b, dist);
    make_perp_line(b, c, a, pbc_b, pbc_c, dist);
    make_perp_line(c, a, b, pca_c, pca_a, dist);
    Vector3f inner_a1, inner_a2;
    Vector3f inner_b1, inner_b2;
    Vector3f inner_c1, inner_c2;
    line_line_intersection(pab_a, pab_b, pca_c, pca_a, inner_a1, inner_a2);
    line_line_intersection(pbc_b, pbc_c, pab_a, pab_b, inner_b1, inner_b2);
    line_line_intersection(pca_c, pca_a, pbc_b, pbc_c, inner_c1, inner_c2);
    const int i = (int)wire_verts.size();
    wire_verts.push_back(a);
    wire_verts.push_back(b);
    wire_verts.push_back(c);
    wire_verts.push_back((inner_a1 + inner_a2)*0.5f);
    wire_verts.push_back((inner_b1 + inner_b2)*0.5f);
    wire_verts.push_back((inner_c1 + inner_c2)*0.5f);
    wire_quads.push_back({i+1, i+4, i+3, i+0});
    wire_quads.push_back({i+2, i+5, i+4, i+1});
    wire_quads.push_back({i+0, i+3, i+5, i+2});
    for (int u = 0; u < 6; ++u) { wire_uvs.push_back(uv); }
  }
  std::ofstream fout(fname);
  if (!fout.is_open()) { return false; }
  for (const Vector3f& vert : wire_verts) {
    fout << "v " << std::setprecision(12) << (double)vert.x() << " " << (double)vert.y() << " " << (double)vert.z() << "\n";
  }
  for (const int& u : wire_uvs) {
    const float uv = (float)((u * 23) % face_check.size());
    fout << "vt " << ((uv + 0.5f) / (float)face_check.size()) << " " << 0.5f << "\n";
  }
  for (Face poly : wire_quads) {
    fout << "f";
    for (int p : poly) {
      fout << " " << (p+1) << "/" << (p+1);
    }
    fout << "\n";
  }
  return true;
}

bool export_cutout(const char* fname, const Verts3D& v3ds, const Planes& planes, float width) {
  std::cout << "Exporting cutout: " << fname << "..." << std::endl;
  //Convert to 2D faces
  std::vector<Verts2D> faces(planes.size());
  float cur_x = 0.0f;
  float cur_y = 0.0f;
  float max_x = 0.0f;
  float max_y = 0.0f;
  for (size_t i = 0; i < planes.size(); ++i) {
    //Project points
    make_2d_projection(v3ds, g_polys[i], planes[i], faces[i]);

    //Figure out a bounding box
    Vector2f minCoord(1e9f, 1e9f);
    Vector2f maxCoord(-1e9f, -1e9f);
    for (Vector2f& v : faces[i]) {
      minCoord = minCoord.cwiseMin(v);
      maxCoord = maxCoord.cwiseMax(v);
    }

    //Transform coordinates
    for (Vector2f& v : faces[i]) {
      v -= minCoord;
      v.x() += cur_x;
      v.y() += cur_y;
    }

    //Advance height to next slot
    cur_y += maxCoord.y() - minCoord.y();
    max_x = std::max(max_x, maxCoord.x() - minCoord.x());

    if (i % 3 == 2) {
      cur_x += max_x;
      max_y = std::max(max_y, cur_y);
      cur_y = 0.0f;
      max_x = 0.0f;
    }
  }

  //Compute the scale factor
  const float scale = width / cur_x;

  std::ofstream fout(fname);
  if (!fout.is_open()) { return false; }
  for (const Verts2D& face : faces) {
    for (const Vector2f& p : face) {
      fout << "v 0 " << double(p.x() * scale) << " " << double(p.y() * scale) << "\n";
    }
  }
  int fIx = 0;
  for (const Verts2D& face : faces) {
    fout << "f";
    for (const Vector2f& p : face) {
      fIx += 1;
      fout << " " << fIx;
    }
    fout << "\n";
  }
  return true;
}

void line_line_intersection(const Vector3f& a1, const Vector3f& a2, const Vector3f& b1, const Vector3f& b2, Vector3f& pa, Vector3f& pb) {
    const Vector3f p13 = a1 - b1;
    const Vector3f p43 = b2 - b1;
    const Vector3f p21 = a2 - a1;
    const float d1343 = p13.dot(p43);
    const float d4321 = p43.dot(p21);
    const float d1321 = p13.dot(p21);
    const float d4343 = p43.dot(p43);
    const float d2121 = p21.dot(p21);
    const float denom = d2121 * d4343 - d4321 * d4321;
    const float numer = d1343 * d4321 - d1321 * d4343;
    const float mua = numer / denom;
    const float mub = (d1343 + d4321 * mua) / d4343;
    pa = a1 + std::min(std::max(mua, 0.0f), 1.0f) * p21;
    pb = b1 + std::min(std::max(mub, 0.0f), 1.0f) * p43;
}
float line_line_dist_sq(const Vector3f& a1, const Vector3f& a2, const Vector3f& b1, const Vector3f& b2) {
    Vector3f pa, pb;
    line_line_intersection(a1, a2, b1, b2, pa, pb);
    return (pa - pb).squaredNorm();
}

float point_triangle_dist_sq(Vector3f p, Vector3f a, Vector3f b, Vector3f c) {
    a -= c;
    b -= c;
    p -= c;

    const float aa = a.dot(a);
    const float ab = a.dot(b);
    const float bb = b.dot(b);
    const float av = -a.dot(p);
    const float bv = -b.dot(p);

    float det = aa * bb - ab * ab;
    float s = ab * bv - bb * av;
    float t = ab * av - aa * bv;

    if (s + t < det) {
        if (s < 0.0f) {
            if (t < 0.0f) {
                if (av < 0.0f) {
                    s = std::clamp(-av / aa, 0.0f, 1.0f);
                    t = 0.0f;
                } else {
                    s = 0.0f;
                    t = std::clamp(-bv / bb, 0.0f, 1.0f);
                }
            } else {
                s = 0.0f;
                t = std::clamp(-bv / bb, 0.0f, 1.0f);
            }
        } else if (t < 0.0f) {
            s = std::clamp(-av / aa, 0.0f, 1.0f);
            t = 0.0f;
        } else {
            const float invDet = 1.0f / det;
            s *= invDet;
            t *= invDet;
        }
    } else {
        if (s < 0.0f) {
            const float tmp0 = ab + av;
            const float tmp1 = bb + bv;
            if (tmp1 > tmp0) {
                const float numer = tmp1 - tmp0;
                const float denom = aa - 2 * ab + bb;
                s = std::clamp(numer / denom, 0.0f, 1.0f);
                t = 1.0f - s;
            } else {
                t = std::clamp(-bv / bb, 0.0f, 1.0f);
                s = 0.0f;
            }
        } else if (t < 0.0f) {
            if (aa + av > ab + bv) {
                const float numer = bb + bv - ab - av;
                const float denom = aa - 2 * ab + bb;
                s = std::clamp(numer / denom, 0.0f, 1.0f);
                t = 1.0f - s;
            } else {
                s = std::clamp(-bv / bb, 0.0f, 1.0f);
                t = 0.0f;
            }
        } else {
            const float numer = bb + bv - ab - av;
            const float denom = aa - 2 * ab + bb;
            s = std::clamp(numer / denom, 0.0f, 1.0f);
            t = 1.0f - s;
        }
    }
    return (a * s + b * t - p).squaredNorm();
}

int petrie_length(const Faces& tris) {
    int sum = 0;
    int a = tris[0][0];
    int b = tris[0][1];
    int c = tris[0][2];
    int start_c = c;
    for (int i = 0; i < 9999; ++i) {
        for (const Face& tri : tris) {
            if (std::find(tri.begin(), tri.end(), a) != tri.end() &&
                std::find(tri.begin(), tri.end(), b) != tri.end() &&
                std::find(tri.begin(), tri.end(), c) == tri.end()) {
                const int d = tri[0] + tri[1] + tri[2] - a - b;
                c = b;
                b = a;
                a = d;
                break;
            }
        }
        sum += 1;
        if (c == start_c) {
            return sum;
        }
    }
    std::cout << "WARNING: Could not find a Petrie polygon." << std::endl;
    return 0;
}

void make_2d_projection(const Verts3D v3ds, const Face& poly, const Plane& plane, Verts2D& v2ds) {
    Matrix3f basis; Vector3f p;
    make_2d_projection(v3ds, poly, plane, v2ds, basis, p);
}

void make_2d_projection(const Verts3D v3ds, const Face& poly, const Plane& plane, Verts2D& v2ds, Matrix3f& basis, Vector3f& p) {
    //Find basis for plane x,y,n
    const Vector3f& n = plane.n;
    const Vector3f x = (v3ds[poly[1]] - v3ds[poly[0]]).normalized();
    const Vector3f y = n.cross(x);
    basis.transpose() << x, y, n;
    p = n * plane.d;

    //Create an array of projected vertices
    const size_t num = poly.size();
    v2ds.resize(num);
    for (size_t i = 0; i < num; ++i) {
        v2ds[i] = (basis * (v3ds[poly[i]] - p)).head<2>();
    }
}

Plane get_plane(const Verts3D& pts, const Face& poly) {
    Vector3d mean = Vector3d::Zero();
    for (int ix : poly) {
        mean += pts[ix].cast<double>();
    }
    mean /= double(poly.size());
    Matrix3d xx = Matrix3d::Zero();
    for (int ix : poly) {
        const Vector3d x = pts[ix].cast<double>() - mean;
        xx += x * x.transpose();
    }
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(xx, Eigen::ComputeFullU);
    const Vector3d n = svd.matrixU().col(2);
    return Plane(n.cast<float>(), (float)n.dot(mean));
}

Vector3f plane_intersection(const Plane& p1, const Plane& p2, const Plane& p3) {
    Matrix3f m; m.transpose() << p1.n, p2.n, p3.n;
    Vector3f d(p1.d, p2.d, p3.d);
    return m.inverse() * d;
}

void y_to_v3ds(const VectorXf& y, Verts3D& verts) {
  verts.clear();
  for (int i = 0; i < y.size(); i += 3) {
    const Eigen::Map<const Vector3f> sub_x(y.data() + i);
    verts.emplace_back(sub_x);
  }
}

void v3ds_to_y(const Verts3D& verts, VectorXf& y) {
  y.resize(verts.size() * 3);
  for (size_t i = 0; i < verts.size(); ++i) {
    const Vector3f& v = verts[i];
    y[i * 3 + 0] = v.x();
    y[i * 3 + 1] = v.y();
    y[i * 3 + 2] = v.z();
  }
}

void x_to_planes(const VectorXf& x, Planes& planes) {
    planes.clear();
    for (int i = 0; i < x.size(); i += 3) {
        const Eigen::Map<const Vector3f> sub_x(x.data() + i);
        planes.emplace_back(sub_x);
    }
}

void planes_to_x(const Planes& planes, VectorXf& x) {
    x.resize(planes.size() * 3);
    for (size_t i = 0; i < planes.size(); ++i) {
        const Vector3f nd = planes[i].n * planes[i].d;
        x[i*3 + 0] = nd.x();
        x[i*3 + 1] = nd.y();
        x[i*3 + 2] = nd.z();
    }
}

void planes_to_v3ds(const Faces& dual_tris, const Planes& planes, Verts3D& verts) {
    verts.clear();
    for (const Face& face : dual_tris) {
        const Plane& plane_a = planes[face[0]];
        const Plane& plane_b = planes[face[1]];
        const Plane& plane_c = planes[face[2]];
        verts.push_back(plane_intersection(plane_a, plane_b, plane_c));
    }
}

void v3ds_to_planes(const Verts3D& pts, const Faces& polys, Planes& planes) {
    planes.clear();
    if (polys[0].size() == 3) {
        for (const Face& poly : polys) {
            planes.emplace_back(pts[poly[0]], pts[poly[1]], pts[poly[2]]);
        }
    } else {
        for (const Face& poly : polys) {
            planes.push_back(get_plane(pts, poly));
        }
    }
}

void x_to_v3ds(const VectorXf& x, const Faces& tris, Verts3D& verts) {
    Planes planes;
    x_to_planes(x, planes);
    planes_to_v3ds(tris, planes, verts);
}

void v3ds_to_x(const Verts3D& pts, const Faces& polys, VectorXf& x) {
    Planes planes;
    v3ds_to_planes(pts, polys, planes);
    planes_to_x(planes, x);
}

void y_to_x(const VectorXf& n, const VectorXf& y, VectorXf& x) {
    x.resize(n.size());
    for (int i = 0; i < y.size(); ++i) {
        const Eigen::Map<const Vector3f> sub_n(n.data() + i*3);
        Eigen::Map<Vector3f> sub_x(x.data() + i*3);
        const float mag = (std::abs(y[i]) > 1e-3f ? y[i] : 1e-3f);
        sub_x = sub_n.normalized() * mag;
    }
}

void trunc_x(VectorXf& x) {
    for (int i = 0; i < x.size(); ++i) {
        x[i] = std::truncf(x[i]);
    }
}

inline float cp_test(const Vector2f& p1, const Vector2f& p2, const Vector2f& p3) {
    return (p2.x() - p1.x()) * (p3.y() - p1.y()) -
           (p2.y() - p1.y()) * (p3.x() - p1.x());
}

bool point_in_polygon(const Vector2f& p, const Verts2D& pts, int& onEdge) {
    static const float epsilon = 1e-6f;
    static const float epsilon2 = 1e-5f;
    onEdge = -1;
    int windingNumber = 0;
    const Vector2f* p2 = &pts[pts.size() - 1];
    for (size_t i = 0; i < pts.size(); i++) {
        const Vector2f* p1 = &pts[i];
        const float d2 = pt_to_line_dist_sq(p, *p1, *p2);
        if (d2 < epsilon && (*p1 - p).squaredNorm() > epsilon2 && (*p2 - p).squaredNorm() > epsilon2) {
            onEdge = (int)i;
            return true;
        }
        if (p1->y() <= p.y()) {
            if (p2->y() > p.y() && cp_test(*p1, *p2, p) > 0.0f) {
                windingNumber++;
            }
        } else {
            if (p2->y() <= p.y() && cp_test(*p1, *p2, p) < 0.0f) {
                windingNumber--;
            }
        }
        p2 = p1;
    }
    return std::abs(windingNumber) == 1;
}

int count_crossings(const Verts3D& v3ds, const Plane& plane, const Face& poly) {
    static const float epsilon = 1e-3f;
    static const float epsilon2 = 1e-10f;
    static Verts2D v2ds;

    //Initialize results
    int crossings = 0;

    //Create an array of projected vertices
    make_2d_projection(v3ds, poly, plane, v2ds);

    //Iterate over all the line segments in order
    const size_t num = poly.size();
    Vector2f a1 = v2ds[0];
    for (size_t i = 1; i < num; ++i) {
        //Get the first line segment
        const Vector2f& a2 = v2ds[i];
        const Vector2f da = a2 - a1;

        Vector2f b1 = v2ds[num - 1];
        for (size_t j = 0; j < i - 1; ++j) {
            //Ignore edge case that should not intersect
            if (i - j == num - 1) {
                b1 = v2ds[j];
                continue;
            }

            //Get the other line segment
            const Vector2f& b2 = v2ds[j];
            const Vector2f db = b1 - b2;

            //Check if the lines intersect
            const Vector2f ba = b1 - a1;
            const float det = da.x()*db.y() - da.y()*db.x();
            if (det * det < epsilon2) {
                //Lines are parallel. Count as crossing if the distance between is too small.
                const float d = (ba - da * (ba.dot(da) / da.squaredNorm())).squaredNorm();
                if (d < epsilon) {
                    crossings += 1;
                }
            } else {
                //Lines are not parallel. Look for the intersection point if it exists.
                const float t = (db.y()*ba.x() - db.x()*ba.y()) / det;
                if (t > -epsilon && t < 1.0f + epsilon) {
                    const float s = (da.x()*ba.y() - da.y()*ba.x()) / det;
                    if (s > -epsilon && s < 1.0f + epsilon) {
                        crossings += 1;
                    }
                }
            }
            b1 = b2;
        }
        a1 = a2;
    }
    return crossings;
}

int count_crossings(const Verts3D& v3ds, const Planes& planes) {
    int crossings = 0;
    for (size_t i = 0; i < g_polys.size(); ++i) {
        crossings += count_crossings(v3ds, planes[i], g_polys[i]);
    }
    return crossings;
}

int count_intersections(const Verts3D& v3ds, const Planes& planes, const Plane& plane, const Face& poly, const Edges& other_edges) {
    static Verts2D v2ds;

    //Initialize results
    int intersections = 0;

    //Create an array of projected vertices
    Matrix3f basis; Vector3f p;
    make_2d_projection(v3ds, poly, plane, v2ds, basis, p);

    //Calculate axis-aligned bounding box for the polygon
    Vector3f aabbMin(FLT_MAX, FLT_MAX, FLT_MAX);
    Vector3f aabbMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    for (int i : poly) {
        aabbMin = aabbMin.cwiseMin(v3ds[i]);
        aabbMax = aabbMax.cwiseMax(v3ds[i]);
    }

    //Iterate over all valid edges
    for (Edge edge : other_edges) {
        //Get edge vertices
        const int ex1 = edge.first;
        const int ex2 = edge.second;
        const Vector3f& p1 = v3ds[ex1];
        const Vector3f& p2 = v3ds[ex2];

        //Check if AABB intersects edge
        if (std::max(p1.x(), p2.x()) < aabbMin.x() ||
            std::max(p1.y(), p2.y()) < aabbMin.y() ||
            std::max(p1.z(), p2.z()) < aabbMin.z() ||
            std::min(p1.x(), p2.x()) > aabbMax.x() ||
            std::min(p1.y(), p2.y()) > aabbMax.y() ||
            std::min(p1.z(), p2.z()) > aabbMax.z()) { continue; }

        //Get plane intersection with edge
        Vector3f in3d;
        if (!plane.intersect(p1, p2, in3d)) continue;
        const Vector2f in2d = (basis * (in3d - p)).head<2>();

        //Check if plane point is inside the polygon
        int edgeIx = -1;
        if (point_in_polygon(in2d, v2ds, edgeIx)) {
            if (edgeIx >= 0 && poly.size() > 3) {
                //If edges are already crossed, don't double count this intersection
                const int px1 = poly[edgeIx];
                const int px2 = poly[(edgeIx + poly.size() - 1) % poly.size()];
                const float lineDistSq = line_line_dist_sq(v3ds[ex1], v3ds[ex2], v3ds[px1], v3ds[px2]);
                if (lineDistSq < 1e-6f) {
                    //The edges intersect, but are they both part of the same crossed face?
                    //Get the plane and compare it with others
                    const Plane intersection_plane(v3ds[ex1], v3ds[ex2], v3ds[px1]);
                    const Vector3f intersection_nd = intersection_plane.n * intersection_plane.d;
                    static const float MIN_DIST_SQ = 1e-8f;
                    float dist = 999.0f;
                    for (const Plane& test_plane : planes) {
                        dist = (intersection_nd - (test_plane.n * test_plane.d)).squaredNorm();
                        if (dist < MIN_DIST_SQ) { break; }
                    }

                    //Matches an existing plane, this must actually be a crossing, not intersection.
                    if (dist < MIN_DIST_SQ) {
                        continue;
                    }
                }
            }
            intersections += 1;
        }
    }
    return intersections;
}

int count_intersections(const Verts3D& v3ds, const Planes& planes) {
    //Create a quick-access set for the polygon
    static std::vector<Edges> poly_sets;
    if (poly_sets.empty()) {
        poly_sets.resize(g_polys.size());
        for (size_t i = 0; i < g_polys.size(); ++i) {
            poly_sets[i].clear();
            const std::vector<int>& poly = g_polys[i];
            std::unordered_set<int> poly_set(poly.begin(), poly.end());
            for (const Edge& edge : g_edges) {
                const int ex1 = edge.first;
                const int ex2 = edge.second;
                if (poly_set.count(ex1) == 0 && poly_set.count(ex2) == 0) {
                    poly_sets[i].push_back(edge);
                }
            }
        }
    }
    int intersections = 0;
    if (g_sym.perms.size() > 0) {
        const std::vector<Permutation>& perms = (g_dual ? g_sym_tri : g_sym).perms;
        for (const Permutation& perm : perms) {
            const int i = perm.rep();
            intersections += perm.order * count_intersections(v3ds, planes, planes[i], g_polys[i], poly_sets[i]);
        }
    } else {
        for (size_t i = 0; i < g_polys.size(); ++i) {
            intersections += count_intersections(v3ds, planes, planes[i], g_polys[i], poly_sets[i]);
        }
    }
    return intersections;
}

float angle_penalty(const Verts3D& verts) {
    float mdp = 0.0f;
    //float mean_dp = 0.0f;
    //float sum = 0.0f;
    for (const Face& poly : g_polys) {
        const size_t n_sides = poly.size();
        for (size_t i = 0; i < n_sides; ++i) {
            const Vector3f p1 = verts[poly[i]];
            const Vector3f p2 = verts[poly[(i+1) % n_sides]];
            const Vector3f p3 = verts[poly[(i+2) % n_sides]];
            const Vector3f d1 = p1 - p2;
            const Vector3f d2 = p3 - p2;
            const float mags = d1.norm() * d2.norm();
            const float dp = std::abs(d1.dot(d2) / mags);
            mdp = std::max(mdp, dp);
            //mean_dp += dp;
            //sum += 1.0f;
        }
    }
    //mean_dp /= sum;
    //return mdp*0.5f + mean_dp*0.5f;
    return mdp;
}

float dist_quality_old(const Verts3D& verts) {
    float min_dist_sq = 99999.0f;
    float max_dist_sq = 0.0f;
    for (size_t i = 0; i < verts.size(); ++i) {
        const Vector3f& p = verts[i];
        for (const Edge& edge : g_edges) {
            if (edge.first == i || edge.second == i) { continue; }
            const Vector3f& a = verts[edge.first];
            const Vector3f& b = verts[edge.second];
            const float dist_sq = pt_to_line_dist_sq(p, a, b);
            min_dist_sq = std::min(min_dist_sq, dist_sq);
            max_dist_sq = std::max(max_dist_sq, dist_sq);
        }
    }
    const float mratio = min_dist_sq / max_dist_sq;
    return std::sqrt(mratio);
}
float dist_quality(const Verts3D& verts) {
    float min_dist_sq = 99999.0f;
    //float sum_log_dist = 0.0f;
    //float num_dist = 0.0f;
    for (size_t i = 0; i < g_edges.size(); ++i) {
        const Edge& e1 = g_edges[i];
        for (size_t j = 0; j < i; ++j) {
            const Edge& e2 = g_edges[j];
            if (e1.first == e2.first || e1.first == e2.second ||
                e1.second == e2.first || e1.second == e2.second) { continue; }
            const float dist_sq = line_line_dist_sq(verts[e1.first], verts[e1.second], verts[e2.first], verts[e2.second]);
            min_dist_sq = std::min(min_dist_sq, dist_sq);
            //if (std::isfinite(dist_sq)) {
            //    sum_log_dist += 0.5f * std::log(dist_sq);
            //    num_dist += 1.0f;
            //}
        } 
    }
    float max_dist_sq = 0.0f;
    for (size_t i = 0; i < verts.size(); ++i) {
        const Vector3f& a = verts[i];
        for (size_t j = 0; j < i; ++j) {
            const Vector3f& b = verts[j];
            max_dist_sq = std::max(max_dist_sq, (a - b).squaredNorm());
        }
    }
    const float mratio = std::sqrt(min_dist_sq / max_dist_sq);
    return mratio;
    //return mratio*0.5f + 0.5f*std::exp(sum_log_dist / num_dist) / std::sqrt(max_dist_sq);
}
float dist_penalty_old(const Verts3D& verts) {
    return 1.0f - dist_quality_old(verts);
}
float dist_penalty(const Verts3D& verts) {
    return 1.0f - dist_quality(verts);
}

float triangle_quality(const Verts3D& verts) {
    if (!g_dual) { return 1.0f; }
    float min_dist_sq = 99999.0f;
    float max_dist_sq = 0.0f;
    for (const Face& tri : g_polys) {
        const int ixA = tri[0];
        const int ixB = tri[1];
        const int ixC = tri[2];
        for (int i = 0; i < (int)verts.size(); ++i) {
            if (i == ixA || i == ixB || i == ixC) { continue; }
            const float dist_sq = point_triangle_dist_sq(verts[i], verts[ixA], verts[ixB], verts[ixC]);
            min_dist_sq = std::min(dist_sq, min_dist_sq);
            max_dist_sq = std::max(dist_sq, max_dist_sq);
        }
    }
    return std::sqrt(min_dist_sq / max_dist_sq);
}
float triangle_penalty(const Verts3D& verts) {
    return 1.0f - triangle_quality(verts);
}

float length_quality(const Verts3D& verts) {
    float min_dist_sq = 99999.0f;
    float max_dist_sq = 0.0f;
    for (const Edge& edge : g_edges) {
        const Vector3f& a = verts[edge.first];
        const Vector3f& b = verts[edge.second];
        const float dist_sq = (a - b).squaredNorm();
        min_dist_sq = std::min(min_dist_sq, dist_sq);
        max_dist_sq = std::max(max_dist_sq, dist_sq);
    }
    return std::sqrt(min_dist_sq / max_dist_sq);
}
float length_penalty(const Verts3D& verts) {
    return 1.0f - length_quality(verts);
}

float plane_penalty(const Planes& planes) {
  //float mdp = 0.0f;
  float pp = 0.0f;
  float sum = 0.0f;
    for (const Face& triLoop : g_tris) {
        const size_t n_sized = triLoop.size();
        for (size_t i = 0; i < n_sized; ++i) {
            const Plane& p1 = planes[triLoop[i]];
            const Plane& p2 = planes[triLoop[(i + 1) % n_sized]];
            //mdp = std::max(mdp, std::abs(p1.n.dot(p2.n)));
            pp += std::abs(p1.n.dot(p2.n));
            sum += 1.0f;
        }
    }
    //return mdp;
    return pp / sum;
}

float q_penalty(const Verts3D& verts, const Planes& planes, int crossings) {
  //return plane_penalty(planes)*0.5f + length_penalty(verts)*0.5f;
  const float dp = (crossings == 0 ? dist_penalty(verts) : dist_penalty_old(verts));
  const float ap = angle_penalty(verts);
  if (g_dual) {
    const float tp = triangle_penalty(verts);
    const float pp = plane_penalty(planes);
    //return std::max(ap, std::max(dp, pp));
    //return std::max(ap, std::max(dp, std::max(pp, tp)));
    return softmax(ap, dp, pp, tp);
  } else {
    const float lp = length_penalty(verts);
    return std::max(ap, std::max(dp, lp));
    //return softmax(ap, dp, lp);
  }
}

float objective_cross_int_qlim(const Verts3D& v3ds, const Planes& planes) {
    if (angle_penalty(v3ds) > 0.999f || dist_penalty(v3ds) > 0.999f) {
        return 99999.0f;
    }
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_crossings * 100 + num_intersections);
}

float objective_cross_int_q(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    const float q = q_penalty(v3ds, planes, num_crossings);
    return float(num_crossings * 100 + num_intersections) + q;
}

float objective_cross_int(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_crossings * 100 + num_intersections);
}

float objective_cross_zint_q(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    if (num_crossings == 0) {
        const int num_intersections = count_intersections(v3ds, planes);
        return float(num_intersections) + q_penalty(v3ds, planes, num_crossings);
    } else {
        return float(num_crossings * 1000) + q_penalty(v3ds, planes, num_crossings);
    }
}

float objective_cross_zint(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    if (num_crossings == 0) {
        const int num_intersections = count_intersections(v3ds, planes);
        return float(num_intersections);
    } else {
        return float(num_crossings * 1000);
    }
}
float objective_cross(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    return float(num_crossings);
}

float objective_int_cross_q(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    const float q = q_penalty(v3ds, planes, num_crossings);
    return float(num_intersections * 100 + num_crossings) + q;
}

float objective_int_cross(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_intersections * 100 + num_crossings);
}

float objective_int_zcross_q(const Verts3D& v3ds, const Planes& planes) {
    const int num_intersections = count_intersections(v3ds, planes);
    if (num_intersections == 0) {
        const int num_crossings = count_crossings(v3ds, planes);
        const float q = q_penalty(v3ds, planes, num_crossings);
        return float(num_crossings) + q;
    } else {
        return float(num_intersections * 1000);
    }
}

float objective_int_zcross(const Verts3D& v3ds, const Planes& planes) {
    const int num_intersections = count_intersections(v3ds, planes);
    if (num_intersections == 0) {
        const int num_crossings = count_crossings(v3ds, planes);
        return float(num_crossings);
    } else {
        return float(num_intersections * 1000);
    }
}

float objective_int_qlim(const Verts3D& v3ds, const Planes& planes) {
    const float q = length_penalty(v3ds);
    //const float q = dist_penalty_old(v3ds);
    //const float q = std::max(q_penalty(v3ds), plane_penalty(planes));
    if (q >= 0.9999f) {
      return 1e12f;
    }
    //const float q = softmax(q_penalty(v3ds), plane_penalty(planes));
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_intersections);
}

float objective_int_q(const Verts3D& v3ds, const Planes& planes) {
    //if (plane_penalty(planes) > 0.9999f) { return 1.0f; }
    //if (angle_penalty(v3ds) > 0.9999f) { return 1.0f; }
    const float q = q_penalty(v3ds, planes, 0);
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_intersections) + q;
}

float objective_int(const Verts3D& v3ds, const Planes& planes) {
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_intersections);
}

float objective_sum_qlim(const Verts3D& v3ds, const Planes& planes) {
    if (angle_penalty(v3ds) > 0.999f || dist_penalty(v3ds) > 0.999f) {
        return 1e12f;
    }
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_intersections + num_crossings);
}

float objective_sum_q(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    const float q = q_penalty(v3ds, planes, num_crossings);
    return float(num_intersections + num_crossings) + q;
}

float objective_sum(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_intersections + num_crossings);
}

float objective_max(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    return float(std::max(num_intersections, num_crossings));
}

float objective_wsum(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_intersections + num_crossings * 2);
}

float objective_wsum2(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    return float(num_intersections * 2 + num_crossings);
}

float objective_wsum_q(const Verts3D& v3ds, const Planes& planes) {
    const int num_crossings = count_crossings(v3ds, planes);
    const int num_intersections = count_intersections(v3ds, planes);
    const float q = q_penalty(v3ds, planes, num_crossings);
    return float(num_intersections + num_crossings * 2) + q;
}
