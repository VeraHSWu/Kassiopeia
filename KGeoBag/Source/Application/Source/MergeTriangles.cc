/*
 * MergeTriangles.hh
 *
 *  Created on: 05 December 2022
 *      Author: Jan Behrens
 */

#include "KGCoreMessage.hh"
#include "KGInterfaceBuilder.hh"
#include "KGTriangle.hh"
#include "KGRectangle.hh"

#include "../../Shapes/External/Source/happly.h"
#include "../../Shapes/External/Source/stl_reader.h"

#include <utility>

#ifdef USE_OPENMP
#include <omp.h>
#endif

using namespace KGeoBag;
using namespace katrin;
using namespace std;

class KGTriangleMerger
{
public:

    /// read STL file into list of KGTriangle's
    void readStlFile(const string& aFilename, size_t max_faces = 0);

    /// read PLY file into list of KGTriangle's and KGRectangle's
    void readPlyFile(const string& aFilename,  size_t max_faces = 0);

    /// write PLY file from list of KGTriangle's and KGRectangles's
    void writePlyFile(const string& aFilename, happly::DataFormat aFormat = happly::DataFormat::Binary);

    /// try to merge KGTriangle pairs into single KGRectangle; optionally check surface normals
    void mergeFaces(bool aCheckNormals = true, bool checkDistance = true);

public:
    KGTriangleMerger() = default;
    ~KGTriangleMerger() = default;

    std::vector<KGTriangle>& GetTriangles() { return fTriangles; }
    std::vector<KGRectangle>& GetRectangles() { return fRectangles; }
    std::vector<std::array<double, 3>>& GetVertices() { return fVertices; }
    std::vector<std::vector<size_t>>& GetFaces() { return fFaces; }

protected:
    /// create KGTriangle from STL mesh element
    template<typename ValueT, typename IndexT>
    KGTriangle GetTriangle(stl_reader::StlMesh<ValueT, IndexT>& mesh, size_t index, double scale = 1.) const;

    /// create KGRectangle from PLY mesh element
    template<typename ValueT, typename IndexT>
    KGRectangle GetRectangle(const std::vector<std::array<ValueT, 3>>& vertices, const std::vector<IndexT>& indices, double scale = 1.) const;

    /// create KGTriangle from PLY mesh element
    template<typename ValueT, typename IndexT>
    KGTriangle GetTriangle(const std::vector<std::array<ValueT, 3>>& vertices, const std::vector<IndexT>& indices, double scale = 1.) const;

    /// show progress bar on terminal
    std::ostream& progressBar(int value, int max_value = 0, int length = 50, std::ostream& = std::cout) const;

    /// get vertex index for given point (avoids duplicate vertices)
    size_t getVertexIndex(const KThreeVector& aVector);

    /// compare scalar equality
    inline bool nearlyEqual(const double& a1, const double& a2, double epsilon = 1.e-10f) const;

    /// compare vector equality
    inline bool nearlyEqual(const KThreeVector& a1, const KThreeVector& a2, double epsilon = 1.e-10f) const;

private:
    std::string fInputFilename;
    std::vector<KGTriangle> fTriangles;
    std::vector<KGRectangle> fRectangles;

    struct ArrayHasher {
        std::size_t operator()(const std::array<double, 3>& a) const {
            std::size_t h = 0;

            for (auto e : a) {
                h ^= std::hash<double>{}(e)  + 0x9e3779b9 + (h << 6) + (h >> 2);
            }
            return h;
        }
    };

    std::unordered_map<std::array<double, 3>, size_t, ArrayHasher> fVertexIndexMap;
    std::vector<std::array<double, 3>> fVertices;
    std::vector<std::vector<size_t>> fFaces;
};

/// create KGTriangle from STL mesh element
template<typename ValueT, typename IndexT>
KGTriangle KGTriangleMerger::GetTriangle(stl_reader::StlMesh<ValueT, IndexT>& mesh, size_t index, double scale) const
{
    KThreeVector p[3], n;
    for (size_t icorner = 0; icorner < 3; ++icorner) {
        const double* c = mesh.tri_corner_coords(index, icorner);
        p[icorner] = c;
    }
    n = mesh.tri_normal(index);

    auto tri = KGTriangle(p[0] * scale, p[1] * scale, p[2] * scale);
    if (n.Dot(tri.GetN3()) < 0) {
        // normals are flipped, so change order of points
        tri = KGTriangle(p[2] * scale, p[1] * scale, p[0] * scale);
    }

    return tri;
}

/// create KGTriangle from PLY mesh element
template<typename ValueT, typename IndexT>
KGTriangle KGTriangleMerger::GetTriangle(const std::vector<std::array<ValueT, 3>>& vertices, const std::vector<IndexT>& indices, double scale) const
{
    KThreeVector p[3], n;
    assert(indices.size() == 3);

    for (auto icorner = 0; icorner < 3; icorner++) {
        p[icorner] = vertices[indices[icorner]];
    }

    auto tri = KGTriangle(p[0] * scale, p[1] * scale, p[2] * scale);

    return tri;
}

/// create KGRectangle from PLY mesh element
template<typename ValueT, typename IndexT>
KGRectangle KGTriangleMerger::GetRectangle(const std::vector<std::array<ValueT, 3>>& vertices, const std::vector<IndexT>& indices, double scale) const
{
    KThreeVector p[4], n;
    assert(indices.size() == 4);

    for (auto icorner = 0; icorner < 4; icorner++) {
        p[icorner] = vertices[indices[icorner]];
    }

    auto rect = KGRectangle(p[0] * scale, p[1] * scale, p[2] * scale, p[3] * scale);

    return rect;
}

/// show progress bar on terminal
std::ostream& KGTriangleMerger::progressBar(int value, int max_value, int length, std::ostream& strm) const
{
    const int progress = max_value > 0 ? length * (float) (value+1) / max_value : length;

    strm << "\r  ";
    for (int j = 0; j < length; j++)
        strm << (j <= progress ? "#" : ".");
    strm << "  [" << 2 * progress << "%]" << std::flush;
    return strm;
}

/// get vertex index for given point (avoids duplicate vertices)
size_t KGTriangleMerger::getVertexIndex(const KThreeVector& aVector)
{
    const auto& point = aVector.AsArray();
    size_t index = 0;

    // find index of already existing vertex
    if (fVertexIndexMap.find(point) != fVertexIndexMap.end()) {
        index = fVertexIndexMap[point];
    }
    // otherwise, append vertex and return new index
    else {
        index = fVertices.size();  // new index at end of list
        fVertices.push_back(point);
        fVertexIndexMap[point] = index;
    }

    return index;
}

/// compare scalar equality
inline bool KGTriangleMerger::nearlyEqual(const double& a1, const double& a2, double epsilon) const
{
    return (fabs(a1 - a2) < epsilon) ? true : false;
}

/// compare vector equality
inline bool KGTriangleMerger::nearlyEqual(const KThreeVector& a1, const KThreeVector& a2, double epsilon) const
{
    return (fabs(a1.GetX() - a2.GetX()) < epsilon &&
            fabs(a1.GetY() - a2.GetY()) < epsilon &&
            fabs(a1.GetZ() - a2.GetZ()) < epsilon) ? true : false;
}

/// read STL file into list of KGTriangle's
void KGTriangleMerger::readStlFile(const string& aFilename, size_t max_faces)
{
    // Adapted from https://github.com/sreiter/stl_reader
    try {
        // read input file
        coremsg_debug("reading elements from STL file <" << aFilename << ">" << eom);
        stl_reader::StlMesh <double, unsigned int> mesh(aFilename.c_str());

        const auto num_tris = mesh.num_tris();

        fTriangles.clear();
        fTriangles.reserve(num_tris);

        // convert all faces into triangles
        for(size_t itri = 0; itri < num_tris; ++itri) {
            if (max_faces > 0 && itri >= max_faces)
                break;

            KGTriangle tri = GetTriangle(mesh, itri);
            fTriangles.emplace_back(tri);
        }
    }
    catch (std::exception &e) {
        coremsg(eError) << "could not read from file <" << aFilename << ">: " << e.what() << eom;
        throw;
    }

    fInputFilename = aFilename;
}

/// read PLY file into list of KGTriangle's and KGRectangle's
void KGTriangleMerger::readPlyFile(const string& aFilename, size_t max_faces)
{
    // Adapted from https://github.com/nmwsharp/happly
    try {
        // read input file
        coremsg_debug("reading elements from PLY file <" << aFilename << ">" << eom);
        happly::PLYData mesh(aFilename.c_str());

        const std::vector<std::array<double, 3>> vertex_positions = mesh.getVertexPositions();
        const std::vector<std::vector<size_t>> face_indices = mesh.getFaceIndices<size_t>();

        const auto num_faces = face_indices.size();

        fTriangles.clear();
        fTriangles.reserve(num_faces);

        fRectangles.clear();
        fRectangles.reserve(num_faces);

        // convert all faces into triangles and rectangles (higher orders are ignored)
        for (size_t iface = 0; iface < num_faces; ++iface) {
            if (max_faces > 0 && iface >= max_faces)
                break;

            const size_t n_dim = face_indices[iface].size();
            switch (n_dim) {
                case 3:
                {
                    auto tri = GetTriangle(vertex_positions, face_indices[iface]);
                    fTriangles.emplace_back(tri);
                    break;
                }
                case 4:
                {
                    auto rect = GetRectangle(vertex_positions, face_indices[iface]);
                    fRectangles.emplace_back(rect);
                    break;
                }
                default:
                    coremsg(eError) << n_dim << "-sided face in ply file <" << aFilename << "> is not supported (index " << iface << ")";
            }
        }
    }
    catch (std::exception &e) {
        coremsg(eError) << "could not read from file <" << aFilename << ">: " << e.what() << eom;
        throw;
    }

    fInputFilename = aFilename;
}

/// write PLY file from list of KGTriangle's and KGRectangles's
void KGTriangleMerger::writePlyFile(const string& aFilename, happly::DataFormat aFormat)
{
    // prepare internal buffers
    fVertexIndexMap.clear();
    fVertices.clear();
    fFaces.clear();

    fFaces.reserve(fTriangles.size() + fRectangles.size());

    // create list of indices for all triangle surfaces
    for (auto & tri : fTriangles) {
        std::vector<size_t> indices;
        indices.push_back(getVertexIndex(tri.GetP0()));
        indices.push_back(getVertexIndex(tri.GetP1()));
        indices.push_back(getVertexIndex(tri.GetP2()));

        fFaces.push_back(indices);
    }

    // create list of indices for all rectangle surfaces
    for (auto & rect : fRectangles) {
        std::vector<size_t> indices;
        indices.push_back(getVertexIndex(rect.GetP0()));
        indices.push_back(getVertexIndex(rect.GetP1()));
        indices.push_back(getVertexIndex(rect.GetP2()));
        indices.push_back(getVertexIndex(rect.GetP3()));

        fFaces.push_back(indices);
    }

    // write output file
    happly::PLYData mesh;
    mesh.addVertexPositions(fVertices);
    mesh.addFaceIndices(fFaces);

    try {
        coremsg_debug("writing elements to PLY file <" << aFilename << ">" << eom);
        mesh.write(aFilename, aFormat);
    }
    catch (std::exception &e) {
        coremsg(eError) << "could not write to file <" << aFilename << ">: " << e.what() << eom;
        throw;
    }
}

/// try to merge KGTriangle pairs into single KGRectangle; optionally check surface normals
void KGTriangleMerger::mergeFaces(bool aCheckNormals, bool checkDistance)
{
    // prepare internal buffers
    std::vector<KGTriangle> triangles_out;
    std::vector<KGRectangle> rectangles_out;
    std::vector<std::pair<KThreeVector,double>> bounding_spheres;  // center + radius
    std::vector<bool> merged_faces;  // flag

    const size_t total_faces = fTriangles.size() + fRectangles.size();
    size_t merged_count = 0;

    triangles_out.reserve(total_faces);
    rectangles_out.reserve(total_faces);
    bounding_spheres.reserve(total_faces);
    merged_faces.reserve(total_faces);

    // fill internal buffers
    std::fill_n(merged_faces.begin(), total_faces, false);

    if (checkDistance) {
        std::transform(fTriangles.cbegin(), fTriangles.cend(), bounding_spheres.begin(),
                       [](const KGTriangle& t)
        {
            KThreeVector C = t.GetP0() + (t.GetA() * t.GetN1() + t.GetB() * t.GetN2()) / 3.;  // centroid
            double R = sqrt(fmax((t.GetP0() - C).MagnitudeSquared(),
                                 fmax((t.GetP1() - C).MagnitudeSquared(),
                                      (t.GetP2() - C).MagnitudeSquared())));
            return std::make_pair(C, R);
        });
    }

    // outer loop over all triangles from input
    for (size_t i = 0; i < fTriangles.size(); ++i) {

        if (i % 100 == 0)
            progressBar(i, total_faces) << " (scanned: " << i << ", merged: " << merged_count << ")" << std::flush;

        if (merged_faces[i])
            continue;

        auto& tri1 = fTriangles[i];
        KThreeVector& p1 = bounding_spheres[i].first;
        double& r1 = bounding_spheres[i].second;

        volatile bool is_merged = false;

#ifdef USE_OPENMP
#pragma omp parallel for shared(is_merged)
#endif
        // inner loop over all REMAINING (not yet compared) triangles
        for (size_t j = i+1; j < fTriangles.size(); ++j) {

            if (is_merged || merged_faces[j])
                continue;

            auto& tri2 = fTriangles[j];
            KThreeVector& p2 = bounding_spheres[j].first;
            double& r2 = bounding_spheres[j].second;

            // check (1): triangles must be close to each other
            if (checkDistance && (p1 - p2).MagnitudeSquared() > (r1 + r2) * (r1 + r2)) {
                continue;
            }

            // check (2): surface normals of triangles must match
            if (aCheckNormals && ! nearlyEqual(tri1.GetN3(), tri2.GetN3())) {
                continue;
            }

            // try to find rectangle that matches all 6 triangle vertices
            std::array<KThreeVector, 4> points;

            //--------------
            // 0-1-2 + 0-1-2
            /**/ if (nearlyEqual(tri1.GetP0(), tri2.GetP2()) && nearlyEqual(tri1.GetP2(), tri2.GetP0()))
                points = { tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP1() };

            // 0-1-2 + 1-2-0
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP0()) && nearlyEqual(tri1.GetP2(), tri2.GetP1()))
                points = { tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP2() };

            // 0-1-2 + 2-0-1
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP1()) && nearlyEqual(tri1.GetP2(), tri2.GetP2()))
                points = { tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP0() };

            // 0-1-2 + 2-1-0
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP0()) && nearlyEqual(tri1.GetP2(), tri2.GetP2()))
                points = { tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP1() };

            // 0-1-2 + 0-2-1
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP1()) && nearlyEqual(tri1.GetP2(), tri2.GetP0()))
                points = { tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP2() };

            // 0-1-2 + 1-0-2
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP2()) && nearlyEqual(tri1.GetP2(), tri2.GetP1()))
                points = { tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP0() };

            //--------------
            // 1-2-0 + 0-1-2
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP2()) && nearlyEqual(tri1.GetP0(), tri2.GetP0()))
                points = { tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP1() };

            // 1-2-0 + 1-2-0
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP0()) && nearlyEqual(tri1.GetP0(), tri2.GetP1()))
                points = { tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP2() };

            // 1-2-0 + 2-0-1
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP1()) && nearlyEqual(tri1.GetP0(), tri2.GetP2()))
                points = { tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP0() };

            // 1-2-0 + 2-1-0
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP0()) && nearlyEqual(tri1.GetP0(), tri2.GetP2()))
                points = { tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP1() };

            // 1-2-0 + 0-2-1
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP1()) && nearlyEqual(tri1.GetP0(), tri2.GetP0()))
                points = { tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP2() };

            // 1-2-0 + 1-0-2
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP2()) && nearlyEqual(tri1.GetP0(), tri2.GetP1()))
                points = { tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP0() };

            //--------------
            // 2-0-1 + 0-1-2
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP2()) && nearlyEqual(tri1.GetP1(), tri2.GetP0()))
                points = { tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP1() };

            // 2-0-1 + 1-2-0
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP0()) && nearlyEqual(tri1.GetP1(), tri2.GetP1()))
                points = { tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP2() };

            // 2-0-1 + 2-0-1
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP1()) && nearlyEqual(tri1.GetP1(), tri2.GetP2()))
                points = { tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP0() };

            // 2-0-1 + 2-1-0
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP0()) && nearlyEqual(tri1.GetP1(), tri2.GetP2()))
                points = { tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP1() };

            // 2-0-1 + 0-2-1
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP1()) && nearlyEqual(tri1.GetP1(), tri2.GetP0()))
                points = { tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP2() };

            // 2-0-1 + 1-0-2
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP2()) && nearlyEqual(tri1.GetP1(), tri2.GetP1()))
                points = { tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP0() };

            //--------------
            // 0-2-1 + 0-1-2
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP2()) && nearlyEqual(tri1.GetP1(), tri2.GetP0()))
                points = { tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP1() };

            // 0-2-1 + 1-2-0
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP0()) && nearlyEqual(tri1.GetP1(), tri2.GetP1()))
                points = { tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP2() };

            // 0-2-1 + 2-0-1
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP1()) && nearlyEqual(tri1.GetP1(), tri2.GetP2()))
                points = { tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP0() };

            // 0-2-1 + 2-1-0
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP0()) && nearlyEqual(tri1.GetP1(), tri2.GetP2()))
                points = { tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP1() };

            // 0-2-1 + 0-2-1
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP1()) && nearlyEqual(tri1.GetP1(), tri2.GetP0()))
                points = { tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP2() };

            // 0-2-1 + 1-0-2
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP2()) && nearlyEqual(tri1.GetP1(), tri2.GetP1()))
                points = { tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP0() };

            //--------------
            // 1-0-2 + 0-1-2
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP2()) && nearlyEqual(tri1.GetP2(), tri2.GetP0()))
                points = { tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP1() };

            // 1-0-2 + 1-2-0
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP0()) && nearlyEqual(tri1.GetP2(), tri2.GetP1()))
                points = { tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP2() };

            // 1-0-2 + 2-0-1
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP1()) && nearlyEqual(tri1.GetP2(), tri2.GetP2()))
                points = { tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP0() };

            // 1-0-2 + 2-1-0
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP0()) && nearlyEqual(tri1.GetP2(), tri2.GetP2()))
                points = { tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP1() };

            // 1-0-2 + 0-2-1
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP1()) && nearlyEqual(tri1.GetP2(), tri2.GetP0()))
                points = { tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP2() };

            // 1-0-2 + 1-0-2
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP2()) && nearlyEqual(tri1.GetP2(), tri2.GetP1()))
                points = { tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP0() };

            //--------------
            // 2-1-0 + 0-1-2
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP2()) && nearlyEqual(tri1.GetP0(), tri2.GetP0()))
                points = { tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP1() };

            // 2-1-0 + 1-2-0
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP0()) && nearlyEqual(tri1.GetP0(), tri2.GetP1()))
                points = { tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP2() };

            // 2-1-0 + 2-0-1
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP1()) && nearlyEqual(tri1.GetP0(), tri2.GetP2()))
                points = { tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP0() };

            // 2-1-0 + 2-1-0
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP0()) && nearlyEqual(tri1.GetP0(), tri2.GetP2()))
                points = { tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP1() };

            // 2-1-0 + 0-2-1
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP1()) && nearlyEqual(tri1.GetP0(), tri2.GetP0()))
                points = { tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP2() };

            // 2-1-0 + 1-0-2
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP2()) && nearlyEqual(tri1.GetP0(), tri2.GetP1()))
                points = { tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP0() };

            //--------------
            else
                continue;

            // check (3): sides of rectangle must be parallel
            if (! (nearlyEqual(points[1] - points[0], points[2] - points[3]) &&
                   nearlyEqual(points[3] - points[0], points[2] - points[1])) ) {
                continue;
            }

            KGRectangle rect;
            try {
                rect = KGRectangle(points[0], points[1], points[2], points[3]);
            }
            catch (...) {
                coremsg(eDebug) << "merged rectangle is not parallel for triangles " << i << "," << j << eom;
                continue;
            }

            // check(4): normal of rectangle must match triangle's
            if (aCheckNormals & ! nearlyEqual(tri1.GetN3(), rect.GetN3())) {
                rect.FlipSurface();  // try inverted normal

                if (! nearlyEqual(tri1.GetN3(), rect.GetN3())) {
                    coremsg(eDebug) << "merged rectangle does not match normals of triangles " << i << "," << j << eom;
                    continue;
                }
            }

#ifdef USE_OPENMP
#pragma omp critical
#endif
            // mark input surfaces as merged -> will be ignored in following loops
            if (! is_merged) {
                merged_faces[i] = true;
                merged_faces[j] = true;
                is_merged = true;

                rectangles_out.push_back(rect);
                merged_count++;
            }
        }

#ifdef USE_OPENMP
#pragma omp critical
#endif
        // if no match was found, simply copy triangle from input
        if (! is_merged) {
            triangles_out.push_back(tri1);
        }
    }

    // rectangles from input are just copied over
    for (size_t i = 0; i < fRectangles.size(); ++i) {

        if (i % 100 == 0) {
            size_t ifull = i + fTriangles.size();
            progressBar(ifull, total_faces) << " (scanned: " << ifull << ", merged: " << merged_count << ")" << std::flush;
        }

        auto& rect = fRectangles[i];

        rectangles_out.push_back(rect);
    }

    progressBar(total_faces) << " (scanned: " << total_faces << ", merged: " << merged_count << ")" << std::endl;  // finalize

    coremsg(eInfo) << merged_count << " triangles were merged" << eom;

    // replace buffers with merge result
    fTriangles = std::move(triangles_out);
    fRectangles = std::move(rectangles_out);
}

/// application
int main(int argc, char** argv)
{
    if (argc < 3) {
        cout << "usage: ./MergeTriangles [-v] [--ascii] [--only-rectangles] [--ignore-normals] [--ignore-distance] <input_file> <output_file> [max_faces]" << endl;
        cout << "  Merge triangles from STL or PLY mesh file into rectangls where possible; save result as PLY file." << endl;
        return -1;
    }

    bool verbose = true;
    bool useAscii = false;
    bool checkNormals = true;
    bool checkDistance = true;
    bool saveOnlyRectangles = false;
    int maxFaces = 0;
    string tInputFilename;
    string tOutputFilename;

    // parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (argv[i] == string("--verbose") || argv[i] == string("-v")) {
            verbose = true;
            continue;
        }
        if (argv[i] == string("--ascii") || argv[i] == string("-a")) {
            useAscii = true;
            continue;
        }
        if (argv[i] == string("--ignore-normals") || argv[i] == string("-N")) {
            checkNormals = false;
            continue;
        }
        if (argv[i] == string("--ignore-distance") || argv[i] == string("-D")) {
            checkDistance = false;
            continue;
        }
        if (argv[i] == string("--only-rectangles") || argv[i] == string("-R")) {
            saveOnlyRectangles = true;
            continue;
        }
        if (tInputFilename.empty()) {
            tInputFilename = argv[i];
            continue;
        }
        if (tOutputFilename.empty()) {
            tOutputFilename = argv[i];
            continue;
        }
        try {
            maxFaces = stoi(argv[i]);
        }
        catch (...) {
            coremsg(eError) << "Unrecognized argument: " << argv[i] << eom;
        }
    }

    KMessageTable::GetInstance().SetTerminalVerbosity(verbose ? KMessageSeverity::eDebugMessage : KMessageSeverity::eInfoMessage);

    // perform conversion
    KGTriangleMerger converter;

    auto tStart = clock();

    coremsg(eNormal) << "Reading input file: " << tInputFilename << eom;
    if (tInputFilename.substr(tInputFilename.size() - 4) == ".stl") {
        converter.readStlFile(tInputFilename, maxFaces);
    }
    else if (tInputFilename.substr(tInputFilename.size() - 4) == ".ply") {
        converter.readPlyFile(tInputFilename, maxFaces);
    }
    else {
        coremsg(eError) << "Unrecognized file format: " << tInputFilename << eom;
    }

    coremsg(eNormal)   << "Input file <" << tInputFilename << "> contains <"
                     << converter.GetTriangles().size() << "> triangles and <"
                     << converter.GetRectangles().size() << "> rectangles." << eom;

    coremsg(eNormal) << "Merging triangles ..."
#ifdef USE_OPENMP
                     << " (OpenMP enabled, " << omp_get_max_threads() << " threads)"
#endif
                     << eom;

    converter.mergeFaces(checkNormals, checkDistance);

    if (saveOnlyRectangles) {
        converter.GetTriangles().clear();
    }

    coremsg(eNormal) << "Writing output file: " << tOutputFilename << " (" << (useAscii ? "ASCII" : "binary") << " format)" << eom;
    converter.writePlyFile(tOutputFilename, useAscii ? happly::DataFormat::ASCII : happly::DataFormat::Binary);

    coremsg(eNormal)   << "Output file <" << tOutputFilename << "> contains <"
                     << converter.GetTriangles().size() << "> triangles and <"
                     << converter.GetRectangles().size() << "> rectangles with <"
                     << converter.GetVertices().size() << "> vertices." << eom;

    auto tStop = clock();
    auto tWalltime = ((double) (tStop - tStart)) / CLOCKS_PER_SEC;  // time in seconds

    coremsg(eInfo) << "The conversion took " << tWalltime << " seconds." << eom;

    return 0;
}
