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
    void mergeFaces(bool aCheckNormals = true);

public:
    KGTriangleMerger() = default;
    ~KGTriangleMerger() = default;

    const std::vector<KGTriangle>& GetTriangles() { return fTriangles; }
    const std::vector<KGRectangle>& GetRectangles() { return fRectangles; }
    const std::vector<std::array<double, 3>>& GetVertices() { return fVertices; }
    const std::vector<std::vector<size_t>>& GetFaces() { return fFaces; }

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

    /// compare vector equality
    inline bool nearlyEqual(const KThreeVector& a1, const KThreeVector& a2, double epsilon = 1.e-12f) const;

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

    KGRectangle rect = KGRectangle(p[0] * scale, p[1] * scale, p[2] * scale, p[3] * scale);

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
void KGTriangleMerger::mergeFaces(bool aCheckNormals)
{
    // prepare internal buffers
    std::unordered_map<size_t, bool> merged_faces;
    std::vector<KGTriangle> triangles_out;
    std::vector<KGRectangle> rectangles_out;

    const size_t total_faces = fTriangles.size() + fRectangles.size();
    size_t mismatching_normals = 0;

    // outer loop over all triangles from input
    for (size_t i = 0; i < fTriangles.size(); ++i) {
        if (i % 100 == 0)
            progressBar(i, total_faces) << " (scanned: " << i << ", merged: " << merged_faces.size() << ")" << std::flush;

        if (merged_faces.find(i) != merged_faces.end())
            continue;

        auto& tri1 = fTriangles[i];
        bool is_merged = false;

        // inner loop over all REMAINING (not yet compared) triangles
        for (size_t j = i+1; j < fTriangles.size(); ++j) {

            if (merged_faces.find(j) != merged_faces.end())
                continue;

            auto& tri2 = fTriangles[j];

            // try to find rectangle that matches all 6 triangle vertices
            auto rect = KGRectangle();

            //--------------
            // 0-1-2 + 0-1-2
            /**/ if (nearlyEqual(tri1.GetP0(), tri2.GetP2()) && nearlyEqual(tri1.GetP2(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP1());

            // 0-1-2 + 1-2-0
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP0()) && nearlyEqual(tri1.GetP2(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP2());

            // 0-1-2 + 2-0-1
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP1()) && nearlyEqual(tri1.GetP2(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP0());

            // 0-1-2 + 2-1-0
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP0()) && nearlyEqual(tri1.GetP2(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP1());

            // 0-1-2 + 0-2-1
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP1()) && nearlyEqual(tri1.GetP2(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP2());

            // 0-1-2 + 1-0-2
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP2()) && nearlyEqual(tri1.GetP2(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP1(), tri1.GetP2(), tri2.GetP0());

            //--------------
            // 1-2-0 + 0-1-2
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP2()) && nearlyEqual(tri1.GetP0(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP1());

            // 1-2-0 + 1-2-0
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP0()) && nearlyEqual(tri1.GetP0(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP2());

            // 1-2-0 + 2-0-1
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP1()) && nearlyEqual(tri1.GetP0(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP0());

            // 1-2-0 + 2-1-0
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP0()) && nearlyEqual(tri1.GetP0(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP1());

            // 1-2-0 + 0-2-1
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP1()) && nearlyEqual(tri1.GetP0(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP2());

            // 1-2-0 + 1-0-2
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP2()) && nearlyEqual(tri1.GetP0(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP2(), tri1.GetP0(), tri2.GetP0());

            //--------------
            // 2-0-1 + 0-1-2
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP2()) && nearlyEqual(tri1.GetP1(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP1());

            // 2-0-1 + 1-2-0
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP0()) && nearlyEqual(tri1.GetP1(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP2());

            // 2-0-1 + 2-0-1
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP1()) && nearlyEqual(tri1.GetP1(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP0());

            // 2-0-1 + 2-1-0
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP0()) && nearlyEqual(tri1.GetP1(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP1());

            // 2-0-1 + 0-2-1
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP1()) && nearlyEqual(tri1.GetP1(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP2());

            // 2-0-1 + 1-0-2
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP2()) && nearlyEqual(tri1.GetP1(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP0(), tri1.GetP1(), tri2.GetP0());

            //--------------
            // 0-2-1 + 0-1-2
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP2()) && nearlyEqual(tri1.GetP1(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP1());

            // 0-2-1 + 1-2-0
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP0()) && nearlyEqual(tri1.GetP1(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP2());

            // 0-2-1 + 2-0-1
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP1()) && nearlyEqual(tri1.GetP1(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP0());

            // 0-2-1 + 2-1-0
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP0()) && nearlyEqual(tri1.GetP1(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP1());

            // 0-2-1 + 0-2-1
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP1()) && nearlyEqual(tri1.GetP1(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP2());

            // 0-2-1 + 1-0-2
            else if (nearlyEqual(tri1.GetP0(), tri2.GetP2()) && nearlyEqual(tri1.GetP1(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP0(), tri1.GetP2(), tri1.GetP1(), tri2.GetP0());

            //--------------
            // 1-0-2 + 0-1-2
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP2()) && nearlyEqual(tri1.GetP2(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP1());

            // 1-0-2 + 1-2-0
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP0()) && nearlyEqual(tri1.GetP2(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP2());

            // 1-0-2 + 2-0-1
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP1()) && nearlyEqual(tri1.GetP2(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP0());

            // 1-0-2 + 2-1-0
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP0()) && nearlyEqual(tri1.GetP2(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP1());

            // 1-0-2 + 0-2-1
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP1()) && nearlyEqual(tri1.GetP2(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP2());

            // 1-0-2 + 1-0-2
            else if (nearlyEqual(tri1.GetP1(), tri2.GetP2()) && nearlyEqual(tri1.GetP2(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP1(), tri1.GetP0(), tri1.GetP2(), tri2.GetP0());

            //--------------
            // 2-1-0 + 0-1-2
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP2()) && nearlyEqual(tri1.GetP0(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP1());

            // 2-1-0 + 1-2-0
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP0()) && nearlyEqual(tri1.GetP0(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP2());

            // 2-1-0 + 2-0-1
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP1()) && nearlyEqual(tri1.GetP0(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP0());

            // 2-1-0 + 2-1-0
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP0()) && nearlyEqual(tri1.GetP0(), tri2.GetP2()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP1());

            // 2-1-0 + 0-2-1
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP1()) && nearlyEqual(tri1.GetP0(), tri2.GetP0()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP2());

            // 2-1-0 + 1-0-2
            else if (nearlyEqual(tri1.GetP2(), tri2.GetP2()) && nearlyEqual(tri1.GetP0(), tri2.GetP1()))
                rect = KGRectangle(tri1.GetP2(), tri1.GetP1(), tri1.GetP0(), tri2.GetP0());

            //--------------
            else
                continue;

            if (rect.GetA() == 0 || rect.GetB() == 0) {
                continue;
            }

            // compare surface normals
            if (aCheckNormals) {
                if (! nearlyEqual(tri1.GetN3(), tri2.GetN3())) {
                    mismatching_normals++;
                    continue;
                }
                if (! nearlyEqual(tri1.GetN3(), rect.GetN3())) {
                    rect.FlipSurface();  // try inverted normal

                    if (! nearlyEqual(tri1.GetN3(), rect.GetN3())) {
                        coremsg(eWarning) << "merged rectangle does not match normals of triangles " << i << "," << j << eom;
                        mismatching_normals++;
                        continue;
                    }
                }
            }

            // mark input surfaces as merged -> will be ignored in following loops
            merged_faces[i] = true;
            merged_faces[j] = true;
            is_merged = true;

            rectangles_out.push_back(rect);
            break;
        }

        // if no match was found, simply copy triangle from input
        if (! is_merged) {
            triangles_out.push_back(tri1);
        }
    }

    // rectangles from input are just copied over
    for (size_t i = 0; i < fRectangles.size(); ++i) {
        if (i % 100 == 0) {
            size_t ifull = i + fTriangles.size();
            progressBar(ifull, total_faces) << " (scanned: " << ifull << ", merged: " << merged_faces.size() << ")" << std::flush;
        }

        auto& rect = fRectangles[i];
        rectangles_out.push_back(rect);
    }

    progressBar(total_faces) << " (merged: " << merged_faces.size() << ")" << std::endl;  // finalize

    coremsg(eInfo) << merged_faces.size() << "triangles were merged (" << mismatching_normals << " mismatching normals)" << eom;

    // replace buffers with merge result
    fTriangles = std::move(triangles_out);
    fRectangles = std::move(rectangles_out);
}

/// application
int main(int argc, char** argv)
{

    if (argc < 3) {
        cout << "usage: ./MergeTriangles  [--ascii] [--ignore-normals] [--no-merge] <input_file> <output_file> [max_faces]" << endl;
        cout << "  Merge triangles from STL or PLY mesh file into rectangls where possible; save result as PLY file." << endl;
        return -1;
    }

    bool useAscii = false;
    bool checkNormals = true;
    bool noMerge = false;
    int maxFaces = 0;
    string tInputFilename;
    string tOutputFilename;

    // parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (argv[i] == string("--ascii") || argv[i] == string("-a")) {
            useAscii = true;
            continue;
        }
        if (argv[i] == string("--ignore-normals") || argv[i] == string("-i")) {
            checkNormals = true;
            continue;
        }
        if (argv[i] == string("--no-merge") || argv[i] == string("-n")) {
            noMerge = true;
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

    // perform conversion
    KGTriangleMerger converter;

    if (tInputFilename.substr(tInputFilename.size() - 4) == ".stl") {
        converter.readStlFile(tInputFilename, maxFaces);
    }
    else if (tInputFilename.substr(tInputFilename.size() - 4) == ".ply") {
        converter.readPlyFile(tInputFilename, maxFaces);
    }
    else {
        coremsg(eError) << "Unrecognized file format: " << tInputFilename << eom;
    }

    coremsg(eNormal) << "input file <" << tInputFilename << "> contains <"
                     << converter.GetTriangles().size() << "> triangles and <"
                     << converter.GetRectangles().size() << "> rectangles" << eom;

    if (! noMerge) {
        converter.mergeFaces(checkNormals);
    }

    converter.writePlyFile(tOutputFilename, useAscii ? happly::DataFormat::ASCII : happly::DataFormat::Binary);

    coremsg(eNormal) << "output file <" << tOutputFilename << "> contains <"
                     << converter.GetTriangles().size() << "> triangles and <"
                     << converter.GetRectangles().size() << "> rectangles with <"
                     << converter.GetVertices().size() << "> vertices" << eom;

    return 0;
}
