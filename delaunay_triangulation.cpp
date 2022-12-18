
#include <fstream>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#define DEBUGx
#define CLION

static
void error(const char *message) {
    printf("%s", message);
    exit(1);
}

using K = CGAL::Simple_cartesian<double>;
using Point = K::Point_2;

// mapping from vertex_id to vertices
static std::vector<Point> id2vertex;

#include <random>

#ifndef DEBUG

// get random number generator
static
auto get_URBG() {
    std::random_device rd;
    std::mt19937 g(rd());
    return g;
}

#endif

// parse input file, randomly shuffle parsed points and fill id2vertex
// also add 4 boundary vertices at beginning
static
bool parse_input(const std::string &file) {
    std::ifstream fin(file);
    if (!fin) {
        return false;
    }
    std::string line;
    std::stringstream buffer;
    int num_pts;
    double x, y, z;
    double x_min, x_max, y_min, y_max;

    if (getline(fin, line)) {
        buffer << line;
        buffer >> num_pts;
        buffer.clear();
    } else {
        goto parse_error;
    }
    if (getline(fin, line)) {
        buffer << line;
        buffer >> x >> y >> z;
        buffer.clear();
        if (z != 0.0) {
            goto parse_error;
        }
        x_min = x_max = x;
        y_min = y_max = y;
        id2vertex.emplace_back(x, y);
    } else {
        goto parse_error;
    }
    for (int i = 1; i < num_pts; i++) {
        if (getline(fin, line)) {
            buffer << line;
            buffer >> x >> y >> z;
            buffer.clear();
            if (z != 0.0) {
                goto parse_error;
            }
            if (x < x_min) {
                x_min = x;
            }
            if (x > x_max) {
                x_max = x;
            }
            if (y < y_min) {
                y_min = y;
            }
            if (y > y_max) {
                y_max = y;
            }
            id2vertex.emplace_back(x, y);
        } else {
            goto parse_error;
        }
    }
    fin.close();

#ifndef DEBUG
    std::shuffle(id2vertex.begin(), id2vertex.end(), get_URBG());
#endif

    x_min -= 1;
    x_max += 1;
    y_min -= 1;
    y_max += 1;
    id2vertex.insert(id2vertex.begin(), {x_max, y_max});
    id2vertex.insert(id2vertex.begin(), {x_max, y_min});
    id2vertex.insert(id2vertex.begin(), {x_min, y_min});
    id2vertex.insert(id2vertex.begin(), {x_min, y_max});

    return true;
    parse_error:
    error("Error while parsing.\n");
    return false;
}

using vertex_id = int;
const int VERTEX_ERROR = -1;
using edge_id = int;
const int EDGE_ERROR = -1;
using facet_id = int;
const int FACET_ERROR = -1;

#define NEXT(i) ((i + 1) % 3)
#define PREV(i) ((i + 2) % 3)

// facet class that stores 3 vertices and corresponding facets
class Triangle {
    std::array<vertex_id, 3> vertex_ids{VERTEX_ERROR, VERTEX_ERROR, VERTEX_ERROR};
    std::array<facet_id, 3> facet_ids{FACET_ERROR, FACET_ERROR, FACET_ERROR};
public:
    Triangle() = default;

    Triangle(std::array<vertex_id, 3> &&vertex_ids, std::array<facet_id, 3> &&facet_ids) {
        this->vertex_ids = vertex_ids;
        this->facet_ids = facet_ids;
    }

    [[nodiscard]] vertex_id vertex_id_at(int i) const {
#ifdef DEBUG
        if (i < 0 || i > 3) { error("index error"); }
#endif
        return this->vertex_ids[i];
    }

    [[nodiscard]] Point &vertex(int i) const {
        vertex_id v = this->vertex_id_at(i);
        return id2vertex.at(v);
    }

    [[nodiscard]] facet_id facet_id_at(int i) const {
#ifdef DEBUG
        if (i < 0 || i > 3) { error("index error"); }
#endif
        return this->facet_ids[i];
    }

    [[nodiscard]] Triangle &facet(int i) const;

    void update_facet_id(int replace, int with) {
        for (int i = 0; i < 3; i++) {
            if (this->facet_id_at(i) == replace) {
                this->facet_ids[i] = with;
                return;
            }
        }
#ifdef DEBUG
        printf("replace error");
#endif
    }

    [[nodiscard]] double angle(int i) const {
        auto s0 = this->vertex(i) - this->vertex(NEXT(i));
        s0 /= std::sqrt(s0.squared_length());
        auto s1 = this->vertex(i) - this->vertex(PREV(i));
        s1 /= std::sqrt(s1.squared_length());
        auto angle = scalar_product(s0, s1);
        return std::acos(angle);
    }

    void update_vertex_indices(int bias) {
        for (auto &v: this->vertex_ids) {
            v += bias;
        }
    }
};

#include <map>

// mapping from facet_id to facet
static std::map<facet_id, Triangle> id2facet;

Triangle &Triangle::facet(int i) const {
    facet_id f = this->facet_id_at(i);
    return id2facet.at(f);
}

using Line = K::Line_2;

// locate v_target in current triangulation
// if v_target is in an interior of a facet, return {facet_id, EDGE_ERROR}
// else if v_target is on an edge of a facet, return {facet_id, edge_id}
static
std::pair<facet_id, edge_id> locate_point(vertex_id v_target, facet_id f_start = FACET_ERROR) {
    if (f_start == FACET_ERROR) {
        f_start = id2facet.begin()->first;
    }
    auto facet = id2facet.at(f_start);
    std::vector<edge_id> on;
    for (int i = 0; i < 3; i++) {
        auto result = Line(facet.vertex(i), facet.vertex(NEXT(i))).oriented_side(id2vertex.at(v_target));
        switch (result) {
            case -1: // outer: step to another facet
                return locate_point(v_target, facet.facet_id_at(i));
            case 0: // on edge: delayed report
                on.emplace_back(i);
            case 1: // inner
                continue;
        }
    }
    // all inner
    if (on.empty()) {
        return {f_start, EDGE_ERROR};
    } else {
        // on an edge
        if (on.size() == 1) {
            return {f_start, on[0]};
        }
        printf("multiple point at same place");
        exit(1);
    }
}


/*
     0      |      0
    / \     |     /|\
   /   \    |    / | \
  /  f  \   |   /  v  \
 1-------2  |  1-------2
 */
// split a facet into 3 facets, according to v_new in the interior
static
std::array<facet_id, 3> split_facet_interior(facet_id f_2split, vertex_id v_new) {
    auto facet_2split = id2facet.at(f_2split);
    // collects information
    std::array<vertex_id, 3> v_indices{};
    std::array<facet_id, 3> f_indices{};
    for (int i = 0; i < 3; i++) {
        v_indices[i] = facet_2split.vertex_id_at(i);
        f_indices[i] = facet_2split.facet_id_at(i);
    }
    auto f_new1 = static_cast<facet_id>(id2facet.size());
    facet_id f_new2 = f_new1 + 1;
    // updates neighbor information
    if (f_indices[1] != FACET_ERROR) {
        id2facet.at(f_indices[1]).update_facet_id(f_2split, f_new1);
    }
    if (f_indices[2] != FACET_ERROR) {
        id2facet.at(f_indices[2]).update_facet_id(f_2split, f_new2);
    }
    // updates current information
    auto facet_new0 = Triangle(
            {v_indices[0], v_indices[1], v_new},
            {f_indices[0], f_new1, f_new2}
    );
    auto facet_new1 = Triangle(
            {v_indices[1], v_indices[2], v_new},
            {f_indices[1], f_new2, f_2split}
    );
    auto facet_new2 = Triangle(
            {v_indices[2], v_indices[0], v_new},
            {f_indices[2], f_2split, f_new1}
    );
    id2facet.at(f_2split) = facet_new0;
    id2facet.insert({f_new1, facet_new1});
    id2facet.insert({f_new2, facet_new2});
    return {f_2split, f_new1, f_new2};
}

/*
     0      |     0
    / \     |    /|\
   /   \    |   / | \
  /  f  \   |  / f|1 \
 1---v---3  | 1---v---3
  \  o  /   |  \ o|2 /
   \   /    |   \ | /
    \ /     |    \|/
     2      |     2
 */
// split 2 facet into 4 facets, according to v_new on the edge
static
std::array<facet_id, 4> split_facet_edge(facet_id f_2split, edge_id e, vertex_id v_new) {
    auto facet_2split = id2facet.at(f_2split);
    // collects information
    std::array<vertex_id, 4> v_indices{};
    std::array<facet_id, 4> f_indices{};
    v_indices[1] = facet_2split.vertex_id_at(e);
    v_indices[3] = facet_2split.vertex_id_at(NEXT(e));
    v_indices[0] = facet_2split.vertex_id_at(PREV(e));
    f_indices[0] = facet_2split.facet_id_at(PREV(e));
    f_indices[3] = facet_2split.facet_id_at(NEXT(e));
    auto o_2split = facet_2split.facet_id_at(e);
#ifdef DEBUG
    if (o_2split == FACET_ERROR) {
        printf("boundary issue");
        exit(1);
    }
#endif
    auto facet_other_2split = id2facet.at(o_2split);
    edge_id e_other;
    for (e_other = 0; facet_other_2split.vertex_id_at(e_other) != v_indices[3]; e_other++);
#ifdef DEBUG
    if (e_other < 0 or e_other > 3) {
        printf("e_other");
        exit(1);
    }
#endif
    v_indices[2] = facet_other_2split.vertex_id_at(PREV(e));
    f_indices[1] = facet_other_2split.facet_id_at(NEXT(e));
    f_indices[2] = facet_other_2split.facet_id_at(PREV(e));
    auto f_new1 = static_cast<facet_id>(id2facet.size());
    facet_id f_new2 = f_new1 + 1;
    // updates neighbor information
    if (f_indices[2] != FACET_ERROR) {
        id2facet.at(f_indices[2]).update_facet_id(o_2split, f_new2);
    }
    if (f_indices[3] != FACET_ERROR) {
        id2facet.at(f_indices[3]).update_facet_id(f_2split, f_new1);
    }
    // updates current information
    auto facet_newf = Triangle(
            {v_indices[0], v_indices[1], v_new},
            {f_indices[0], o_2split, f_new1}
    );
    auto facet_newo = Triangle(
            {v_indices[1], v_indices[2], v_new},
            {f_indices[1], f_new2, f_2split}
    );
    auto facet_new1 = Triangle(
            {v_indices[3], v_indices[0], v_new},
            {f_indices[3], f_2split, f_new2}
    );
    auto facet_new2 = Triangle(
            {v_indices[2], v_indices[3], v_new},
            {f_indices[2], f_new1, o_2split}
    );
    id2facet.at(f_2split) = facet_newf;
    id2facet.at(o_2split) = facet_newo;
    id2facet.insert({f_new1, facet_new1});
    id2facet.insert({f_new2, facet_new2});
    return {f_2split, o_2split, f_new1, f_new2};
}


/*
      3             3
    1---3    |    1---3
    |\ o|    |    |f /|
 0  | \ |  2 | 0  | / |  2
    |f \|    |    |/ o|
    2---0    |    2---0
      1             1
 */
// determine if an edge (0,1) is illegal, given the new vertex 2
// if it is, flip the edge to (2,3) and recurrently call this on (1,3) and (0,3)
static
void legalize_edge(facet_id f, edge_id e) {
    auto facet = id2facet.at(f);
    facet_id f_other = facet.facet_id_at(e);
    if (f_other == FACET_ERROR) {
        return;
    }
    // collects information
    auto facet_other = id2facet.at(f_other);
    std::array<vertex_id, 4> v_indices = {};
    std::array<facet_id, 4> f_indices = {};
    v_indices[0] = facet.vertex_id_at(e);
    v_indices[1] = facet.vertex_id_at(NEXT(e));
    v_indices[2] = facet.vertex_id_at(PREV(e));
    f_indices[0] = facet.facet_id_at(NEXT(e));
    f_indices[1] = facet.facet_id_at(PREV(e));
    edge_id e_other;
    for (e_other = 0; facet_other.vertex_id_at(e_other) != v_indices[0]; e_other++);
#ifdef DEBUG
    if (e_other < 0 or e_other > 3) {
        printf("e_other");
        exit(1);
    }
#endif
    v_indices[3] = facet_other.vertex_id_at(NEXT(e_other));
    f_indices[2] = facet_other.facet_id_at(e_other);
    f_indices[3] = facet_other.facet_id_at(NEXT(e_other));
    // computes legality
    double prev_angle = 10, angle;
    for (int i = 0; i < 3; i++) {
        angle = facet.angle(i);
        if (angle < prev_angle) {
            prev_angle = angle;
        }
        angle = facet_other.angle(i);
        if (angle < prev_angle) {
            prev_angle = angle;
        }
    }
    auto facet_new = Triangle({v_indices[2], v_indices[3], v_indices[1]}, {f_other, f_indices[3], f_indices[0]});
    auto facet_other_new = Triangle({v_indices[3], v_indices[2], v_indices[0]}, {f, f_indices[1], f_indices[2]});
    double new_angle = 10;
    for (int i = 0; i < 3; i++) {
        angle = facet_new.angle(i);
        if (angle < new_angle) {
            new_angle = angle;
        }
        angle = facet_other_new.angle(i);
        if (angle < new_angle) {
            new_angle = angle;
        }
    }
    // if quad not convex, return
    if (CGAL::right_turn(id2vertex.at(v_indices[2]), id2vertex.at(v_indices[0]), id2vertex.at(v_indices[3]))
        || CGAL::right_turn(id2vertex.at(v_indices[3]), id2vertex.at(v_indices[1]), id2vertex.at(v_indices[2]))) {
        return;
    }
    // if illegal
    if (new_angle > prev_angle) {
        // updates neighbor information
        if (f_indices[1] != FACET_ERROR) {
            id2facet.at(f_indices[1]).update_facet_id(f, f_other);
        }
        if (f_indices[3] != FACET_ERROR) {
            id2facet.at(f_indices[3]).update_facet_id(f_other, f);
        }
        // updates current information
        id2facet.at(f) = facet_new;
        id2facet.at(f_other) = facet_other_new;
        // calls recurrently
        legalize_edge(f, 1);
        legalize_edge(f_other, 2);
    }
}

// drop 4 boundary vertices, and all relevant facets
static
void drop_boundary() {
    id2vertex = std::vector<Point>(id2vertex.begin() + 4, id2vertex.end());

    auto itr = id2facet.begin();
    while (itr != id2facet.end()) {
        auto &[f, facet] = *itr;
        bool remove = false;
        for (int i = 0; i < 3; i++) {
            if (facet.vertex_id_at(i) < 4) {
                remove = true;
                break;
            }
        }
        if (remove) {
            for (int i = 0; i < 3; i++) {
                auto f_other = facet.facet_id_at(i);
                if (f_other != FACET_ERROR) {
                    id2facet.at(f_other).update_facet_id(f, FACET_ERROR);
                }
            }
            itr = id2facet.erase(itr);
        } else {
            facet.update_vertex_indices(-4);
            itr++;
        }
    }
}

#ifdef DEBUG

static
void print_mesh() {
    for (const auto &pair: id2facet) {
        auto [f, facet] = pair;
        printf("facet %d: vertices", f);
        for (int k = 0; k < 3; k++) {
            printf(" %d", facet.vertex_id_at(k));
        }
        printf(" facets");
        for (int k = 0; k < 3; k++) {
            printf(" %d", facet.facet_id_at(k));
        }
        printf("\n");
    }
}

#endif

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO.h>

using Mesh = CGAL::Surface_mesh<K::Point_3>;

// save mesh in .off format
static
void save_mesh(const std::string &file) {
    Mesh mesh;
    std::vector<Mesh::Vertex_index> vertex_handles;
    for (const auto &p: id2vertex) {
        auto vertex_handle = mesh.add_vertex({p.x(), p.y(), 0.0});
        vertex_handles.emplace_back(vertex_handle);
    }
    for (const auto &pair: id2facet) {
        mesh.add_face(
                vertex_handles[pair.second.vertex_id_at(0)],
                vertex_handles[pair.second.vertex_id_at(1)],
                vertex_handles[pair.second.vertex_id_at(2)]
        );
    }
    if (CGAL::write_mesh(mesh, file)) {
        printf("Result saved at %s\n", file.c_str());
    } else {
        error((std::string("Cannot write to file ") + file + "\n").c_str());
    }
}

#include <filesystem>

namespace fs = std::filesystem;

void delaunay_triangulation(const std::string &in_file, const std::string &out_file) {
    id2vertex.clear();
    id2facet.clear();

    printf("Parsing %s\n", in_file.c_str());
    if (!fs::exists(in_file)) {
        error("File doesn't exist.\n");
    }
    if (!parse_input(in_file)) {
        error("Parsing failed.\n");
    }

#ifdef DEBUG
    printf("The input point set contains %ld points.\n", id2vertex.size());
    for (const auto &point: id2vertex) {
        printf("%lf %lf\n", point.x(), point.y());
    }
#endif

    Triangle init0({0, 1, 2}, {-1, -1, 1});
    Triangle init1({2, 3, 0}, {-1, -1, 0});
    id2facet.insert({0, init0});
    id2facet.insert({1, init1});

    for (int i = 4; i < id2vertex.size(); i++) {
        auto [f_2split, edge] = locate_point(i);
        if (edge == EDGE_ERROR) {
            auto f_indices = split_facet_interior(f_2split, i);
            for (const auto &f: f_indices) {
                legalize_edge(f, 0);
            }
        } else {
#ifdef DEBUG
            printf("__________ edge encountered_________");
#endif
            auto f_indices = split_facet_edge(f_2split, edge, i);
            for (const auto &f: f_indices) {
                legalize_edge(f, 0);
            }
        }
#ifdef DEBUG
        printf("_____after introducint point %d_____\n", i);
        print_mesh();
//        save_mesh(std::to_string(i) + ".off");
#endif
    }
    drop_boundary();
#ifdef DEBUG
    printf("_____after dropping boundary_____\n");
    print_mesh();
#endif
    save_mesh(out_file);
}


int main() {
#ifdef CLION
    fs::path work_path = fs::current_path().parent_path();
    fs::path in_file = work_path / "data/input10000.xyz";
    fs::path out_file = work_path / "output/output10000.off";
#else
    fs::path work_path = fs::current_path();
    fs::path in_file = work_path / "input.xyz";
    fs::path out_file = work_path / "output.off";
#endif
    delaunay_triangulation(in_file, out_file);
}
