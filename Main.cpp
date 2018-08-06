#include <string>
#include <fstream>
#include <iostream>

#include "glm/glm.hpp"



using glm::vec2;
using glm::vec3;
using glm::vec4;



static constexpr bool k_frontHeavy(true); // Shifts the concentration of vertices toward the front of the airfoil
static constexpr bool k_zForward(true); // Rotates airfoil so "forward" is +z and "up" is +y

static int f_naca[4];
static int f_res;
static std::string f_outFilePath;
static std::unique_ptr<double[]> f_symPoints;
static std::unique_ptr<vec3[]> f_vertLocs;
static std::unique_ptr<vec3[]> f_vertNorms;
static std::unique_ptr<int[]> f_indices;



static void genSymetricPoints() {
    f_symPoints.reset(new double[f_res]);

    f_symPoints[0] = 0.0;
    f_symPoints[f_res - 1] = 0.0;
    double xFactor(1.0f / (f_res - 1));
    double t((f_naca[2] * 10 + f_naca[3]) / 100.0);
    for (int i(1); i < f_res - 1; ++i) {
        double x(i * xFactor);
        if constexpr (k_frontHeavy) x = x * x;
        double x2(x * x);
        f_symPoints[i] = 5.0 * t * (
             0.2969 * std::sqrt(x) +
            -0.1260 * x +
            -0.3516 * x2 +
             0.2843 * x * x2 +
            -0.1036 * x2 * x2
        );
    }
}

static vec3 detNorm(const vec3 & p1, const vec3 & p2, const vec3 & p3) {
    vec2 d1(p2.x - p1.x, p2.y - p1.y);
    vec2 d2(p3.x - p2.x, p3.y - p2.y);
    vec2 tangent(glm::normalize(d1 + d2));
    return vec3(-tangent.y, tangent.x, 0.0f);
}

static vec3 detNorm(const vec3 & p1, const vec3 & p2) {
    vec2 d1(p2.x - p1.x, p2.y - p1.y);
    vec2 tangent(glm::normalize(d1));
    return vec3(-tangent.y, tangent.x, 0.0f);
}

static void genLocs() {
    f_vertLocs.reset(new vec3[4 * f_res]);

    float xFactor(1.0f / (f_res - 1));
    f_vertLocs[0] = f_vertLocs[f_res] = vec3(0.0f, 0.0f, -1.0f);
    f_vertLocs[f_res - 1] = f_vertLocs[2 * f_res - 1] = vec3(1.0f, 0.0f, -1.0f);
    for (int i(1); i < f_res - 1; ++i) {
        int ti(i), bi(f_res + i);
        float x(i * xFactor);
        if constexpr (k_frontHeavy) x = x * x;
        f_vertLocs[ti] = vec3(x,  f_symPoints[i], -1.0f);
        f_vertLocs[bi] = vec3(x, -f_symPoints[i], -1.0f);
    }

    for (int i(0), j(2 * f_res); i < 2 * f_res; ++i, ++j) {
        f_vertLocs[j].x = f_vertLocs[i].x;
        f_vertLocs[j].y = f_vertLocs[i].y;
        f_vertLocs[j].z = 1.0f;
    }
}

static void genNorms() {
    f_vertNorms.reset(new vec3[4 * f_res]);

    f_vertNorms[0] = f_vertNorms[f_res] = detNorm(f_vertLocs[f_res + 1], f_vertLocs[0], f_vertLocs[1]);
    f_vertNorms[f_res - 1] = detNorm(f_vertLocs[f_res - 2], f_vertLocs[f_res - 1]);
    f_vertNorms[2 * f_res - 1] = detNorm(f_vertLocs[2 * f_res - 1], f_vertLocs[2 * f_res - 2]);
    for (int i(1); i < f_res - 1; ++i) {
        int ti(i), bi(f_res + i);
        f_vertNorms[ti] = detNorm(f_vertLocs[ti - 1], f_vertLocs[ti], f_vertLocs[ti + 1]);
        f_vertNorms[bi] = detNorm(f_vertLocs[bi + 1], f_vertLocs[bi], f_vertLocs[bi - 1]);
    }

    for (int i(0), j(2 * f_res); i < 2 * f_res; ++i, ++j) {
        f_vertNorms[j] = f_vertNorms[i];
    }
}

static void genIndices() {
    f_indices.reset(new int[3 * 4 * (f_res - 1)]);

    int ii(0);
    for (int i(0); i < f_res - 1; ++i) {
        int li(i);
        int ri(2 * f_res + i);

        f_indices[ii++] = li;
        f_indices[ii++] = ri;
        f_indices[ii++] = ri + 1;

        f_indices[ii++] = ri + 1;
        f_indices[ii++] = li + 1;
        f_indices[ii++] = li;
    }
    for (int i(0); i < f_res - 1; ++i) {
        int li(f_res + i);
        int ri(3 * f_res + i);

        f_indices[ii++] = ri;
        f_indices[ii++] = li;
        f_indices[ii++] = li + 1;

        f_indices[ii++] = li + 1;
        f_indices[ii++] = ri + 1;
        f_indices[ii++] = ri;
    }
}

static void rotate() {
    for (int i(0); i < f_res * 4; ++i) {
        vec3 & loc(f_vertLocs[i]);
        loc = vec3(loc.z, loc.y, -loc.x);

        vec3 & norm(f_vertNorms[i]);
        norm = vec3(norm.z, norm.y, -norm.x);
    }
}

static void writeObj() {
    std::ofstream ofs(f_outFilePath);
    if (!ofs.good()) {
        std::cerr << "Failed to open output file" << std::endl;
        return;
    }
    ofs.precision(std::numeric_limits<float>::max_digits10);

    for (int i(0); i < f_res * 4; ++i) {
        ofs << "v " << f_vertLocs[i].x << " " << f_vertLocs[i].y << " " << f_vertLocs[i].z << std::endl;
    }
    for (int i(0); i < f_res * 4; ++i) {
        ofs << "vn " << f_vertNorms[i].x << " " << f_vertNorms[i].y << " " << f_vertNorms[i].z << std::endl;
    }
    int nIndices(3 * 4 * (f_res - 1));
    for (int i(0); i < nIndices; i += 3) {
        ofs << "f ";
        ofs << f_indices[i + 0] + 1 << "//" << f_indices[i + 0] + 1 << " ";
        ofs << f_indices[i + 1] + 1 << "//" << f_indices[i + 1] + 1 << " ";
        ofs << f_indices[i + 2] + 1 << "//" << f_indices[i + 2] + 1 << std::endl;
    }
}

static void printUsage() {
    std::cout << "Usage: amt <4 digit NACA> <x resolution> <out file path>" << std::endl;
}

static bool parseArgs(int argc, char ** argv) {
    if (argc < 4) {
        printUsage();
        return false;
    }

    std::string naca(argv[1]);
    if (naca.length() != 4) {
        std::cerr << "Invalid NACA" << std::endl;
        return false;
    }
    for (int i(0); i < 4; ++i) {
        if (naca[i] < '0' || naca[i] > '9') {
            std::cerr << "Invalid NACA" << std::endl;
            return false;
        }
        f_naca[i] = naca[i] - '0';
    }

    try {
        f_res = std::stoi(std::string(argv[2]));
    }
    catch (const std::exception &) {
        std::cerr << "Invalid x resolution" << std::endl;
        return false;
    }
    if (f_res < 3) {
        std::cerr << "x resolution must be at least 3" << std::endl;
        return false;
    }
    else if (f_res > 1000) {
        std::cerr << "x resolution may not be greater than 1000" << std::endl;
        return false;
    }

    f_outFilePath = argv[3];

    return true;
}



int main(int argc, char ** argv) {
    if (!parseArgs(argc, argv)) {
        std::exit(-1);
    }

    if (f_naca[0] != 0 || f_naca[1] != 0) {
        std::cerr << "Only symetric NACA are supported at this time" << std::endl;
        std::exit(-1);
    }

    genSymetricPoints();
    genLocs();
    genNorms();
    genIndices();

    if constexpr (k_zForward) rotate();

    writeObj();

    return 0;
}