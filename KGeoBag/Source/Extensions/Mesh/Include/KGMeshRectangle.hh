#ifndef KGMESHRECTANGLE_DEF
#define KGMESHRECTANGLE_DEF

#include "KGMeshElement.hh"
#include "KGRectangle.hh"
#include "KThreeVector.hh"

namespace KGeoBag
{
class KGMeshRectangle : public KGMeshElement
{
  public:
    KGMeshRectangle(const double& a, const double& b, const KGeoBag::KThreeVector& p0, const KGeoBag::KThreeVector& n1,
                    const KGeoBag::KThreeVector& n2);
    KGMeshRectangle(const KGeoBag::KThreeVector& p0, const KGeoBag::KThreeVector& p1,
                    const KGeoBag::KThreeVector& /*p2*/, const KGeoBag::KThreeVector& p3);
    KGMeshRectangle(const KGRectangle& t);
    KGMeshRectangle(const KGMeshRectangle& r);
    ~KGMeshRectangle() override;

    static std::string Name()
    {
        return "mesh_rectangle";
    }

    double Area() const override;
    double Aspect() const override;
    void Transform(const KTransformation& transform) override;

    double NearestDistance(const KGeoBag::KThreeVector& aPoint) const override;
    KGeoBag::KThreeVector NearestPoint(const KGeoBag::KThreeVector& aPoint) const override;
    KGeoBag::KThreeVector NearestNormal(const KGeoBag::KThreeVector& aPoint) const override;
    bool NearestIntersection(const KGeoBag::KThreeVector& aStart, const KGeoBag::KThreeVector& anEnd,
                             KGeoBag::KThreeVector& anIntersection) const override;

    KGPointCloud<KGMESH_DIM> GetPointCloud() const override;
    unsigned int GetNumberOfEdges() const override
    {
        return 4;
    };
    void GetEdge(KThreeVector& start, KGeoBag::KThreeVector& end, unsigned int index) const override;

    //assignment
    inline KGMeshRectangle& operator=(const KGMeshRectangle& r)
    {
        if (&r != this) {
            fA = r.fA;
            fB = r.fB;
            fP0 = r.fP0;
            fN1 = r.fN1;
            fN2 = r.fN2;
        }
        return *this;
    }


    double GetA() const
    {
        return fA;
    }
    double GetB() const
    {
        return fB;
    }
    const KGeoBag::KThreeVector& GetP0() const
    {
        return fP0;
    }
    const KGeoBag::KThreeVector& GetN1() const
    {
        return fN1;
    }
    const KGeoBag::KThreeVector& GetN2() const
    {
        return fN2;
    }
    const KGeoBag::KThreeVector GetN3() const
    {
        return fN1.Cross(fN2);
    }
    const KGeoBag::KThreeVector GetP1() const
    {
        return fP0 + fN1 * fA;
    }
    const KGeoBag::KThreeVector GetP2() const
    {
        return fP0 + fN1 * fA + fN2 * fB;
    }
    const KGeoBag::KThreeVector GetP3() const
    {
        return fP0 + fN2 * fB;
    }
    void GetN1(KThreeVector& n1) const
    {
        n1 = fN1;
    }
    void GetN2(KThreeVector& n2) const
    {
        n2 = fN2;
    }
    void GetN3(KThreeVector& n3) const
    {
        n3 = fN1.Cross(fN2);
    }
    void GetP0(KThreeVector& p0) const
    {
        p0 = fP0;
    }
    void GetP1(KThreeVector& p1) const
    {
        p1 = fP0 + fN1 * fA;
    }
    void GetP2(KThreeVector& p2) const
    {
        p2 = fP0 + fN1 * fA + fN2 * fA;
    }
    void GetP3(KThreeVector& p3) const
    {
        p3 = fP0 + fN2 * fB;
    }

  protected:
    static bool SameSide(const KGeoBag::KThreeVector& point, const KGeoBag::KThreeVector& A,
                         const KGeoBag::KThreeVector& B, const KGeoBag::KThreeVector& C);

    double fA;
    double fB;
    KGeoBag::KThreeVector fP0;
    KGeoBag::KThreeVector fN1;
    KGeoBag::KThreeVector fN2;
};
}  // namespace KGeoBag

#endif /* KGMESHRECTANGLE_DEF */
