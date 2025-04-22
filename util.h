#pragma once
#include "types.h"
#include "plane.h"
#include "symmetry.h"
#include <unordered_set>
#include <vector>
#include <map>
#include <iostream>

std::vector<std::string> split(const std::string& s, char delim);
void ltrim(std::string& s);
void rtrim(std::string& s);
bool ends_with(const std::string& value, const std::string& ending);
int triangle_id(int s1, int s2, int s3);
int triangle_id(const Face& tri);
bool is_finite(const Verts3D& verts);
bool is_finite(const Planes& planes);
void print_faces(const Faces& faces);
bool open_face_file(const char* fname, Faces& tris, bool zero_indexed=false);
bool open_cyclic_file(const char* fname, Faces& tris, bool zero_indexed=false);
bool open_topology(const char* fname, Faces& tris, int ix);
bool verify_topology(const Faces& tris);
void dual_graph(const Faces& faces, Faces& dual_faces, Edges& dual_edges);
void dual_verts(const Faces& faces, const Verts3D& verts, Verts3D& dual_verts);
void save_sample(const char* name, const Planes& planes, const Verts3D& verts, int iter, bool can_save);
void save_dot_graph(const char* fname, const Edges& edges);
void make_edges(const Faces& faces, Edges& edges);
void fix_face_ordering(Faces& polys, const Edges& edges);
bool test_face_ordering(const Faces& polys);
bool import_obj(const char* fname, Verts3D& verts, Faces& polys);
bool export_obj(const char* fname, const Verts3D& verts, const Faces& polys);
bool export_colored_obj(const char* fname, const Verts3D& verts, const Faces& polys);
bool export_wireframe_obj(const char* fname, const Verts3D& verts, const Faces& polys, float extrude);
bool export_cutout(const char* fname, const Verts3D& v3ds, const Planes& planes, float width);
void line_line_intersection(const Vector3f& a1, const Vector3f& a2, const Vector3f& b1, const Vector3f& b2, Vector3f& pa, Vector3f& pb);
float point_triangle_dist_sq(Vector3f p, Vector3f a, Vector3f b, Vector3f c);
int petrie_length(const Faces& tris);

Plane get_plane(const Verts3D& pts, const Face& poly);
void make_2d_projection(const Verts3D v3ds, const Face& poly, const Plane& plane, Verts2D& v2ds);
void make_2d_projection(const Verts3D v3ds, const Face& poly, const Plane& plane, Verts2D& v2ds, Matrix3f& basis, Vector3f& p);
Vector3f plane_intersection(const Plane& p1, const Plane& p2, const Plane& p3);

void y_to_v3ds(const VectorXf& y, Verts3D& verts);
void v3ds_to_y(const Verts3D& verts, VectorXf& y);
void x_to_planes(const VectorXf& x, Planes& planes);
void planes_to_x(const Planes& planes, VectorXf& x);
void v3ds_to_planes(const Verts3D& pts, const Faces& polys, Planes& planes);
void planes_to_v3ds(const Faces& dual_tris, const Planes& planes, Verts3D& verts);
void x_to_v3ds(const VectorXf& x, const Faces& tris, Verts3D& verts);
void v3ds_to_x(const Verts3D& pts, const Faces& polys, VectorXf& x);
void y_to_x(const VectorXf& n, const VectorXf& y, VectorXf& x);
void trunc_x(VectorXf& x);

bool point_in_polygon(const Vector2f& p, const Verts2D& pts, int& onEdge);
int count_crossings(const Verts3D& v3ds, const Plane& plane, const Face& poly);
int count_crossings(const Verts3D& v3ds, const Planes& planes);
int count_intersections(const Verts3D& v3ds, const Planes& planes, const Plane& plane, const Face& poly, const Edges& other_edges);
int count_intersections(const Verts3D& v3ds, const Planes& planes);
float angle_penalty(const Verts3D& verts);
float dist_penalty(const Verts3D& verts);
float length_penalty(const Verts3D& verts);
float plane_penalty(const Planes& planes);
float triangle_penalty(const Verts3D& verts);
float q_penalty(const Verts3D& verts, const Planes& planes, int corssings);

float objective_cross_int_qlim(const Verts3D& v3ds, const Planes& planes);
float objective_cross_int_q(const Verts3D& v3ds, const Planes& planes);
float objective_cross_int(const Verts3D& v3ds, const Planes& planes);
float objective_cross_zint_q(const Verts3D& v3ds, const Planes& planes);
float objective_cross_zint(const Verts3D& v3ds, const Planes& planes);
float objective_cross(const Verts3D& v3ds, const Planes& planes);
float objective_int_cross_q(const Verts3D& v3ds, const Planes& planes);
float objective_int_cross(const Verts3D& v3ds, const Planes& planes);
float objective_int_zcross_q(const Verts3D& v3ds, const Planes& planes);
float objective_int_zcross(const Verts3D& v3ds, const Planes& planes);
float objective_int_qlim(const Verts3D& v3ds, const Planes& planes);
float objective_int_q(const Verts3D& v3ds, const Planes& planes);
float objective_int(const Verts3D& v3ds, const Planes& planes);
float objective_max(const Verts3D& v3ds, const Planes& planes);
float objective_sum_qlim(const Verts3D& v3ds, const Planes& planes);
float objective_sum_q(const Verts3D& v3ds, const Planes& planes);
float objective_sum(const Verts3D& v3ds, const Planes& planes);
float objective_wsum(const Verts3D& v3ds, const Planes& planes);
float objective_wsum2(const Verts3D& v3ds, const Planes& planes);
float objective_wsum_q(const Verts3D& v3ds, const Planes& planes);
