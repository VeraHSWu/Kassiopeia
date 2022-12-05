#ifndef KGTRIANGLE_H_
#define KGTRIANGLE_H_

#include "KGArea.hh"

#include "KThreeVector.hh"

namespace KGeoBag
{
class KGTriangle : public KGArea
{
  public:
    class Visitor
    {
      public:
        Visitor() = default;
        virtual ~Visitor() = default;

        virtual void Visit(KGTriangle*) = 0;
    };

    KGTriangle();
    KGTriangle(const double& a, const double& b, const katrin::KThreeVector& p0, const katrin::KThreeVector& n1, const katrin::KThreeVector& n2);
    KGTriangle(const katrin::KThreeVector& p0, const katrin::KThreeVector& p1, const katrin::KThreeVector& p2);
    KGTriangle(const KGTriangle&);
    KGTriangle& operator=(const KGTriangle&);

    ~KGTriangle() override = default;

    void AreaInitialize() const override {}
    void AreaAccept(KGVisitor* aVisitor) override;
    bool AreaAbove(const katrin::KThreeVector& aPoint) const override;
    katrin::KThreeVector AreaPoint(const katrin::KThreeVector& aPoint) const override;
    katrin::KThreeVector AreaNormal(const katrin::KThreeVector& aPoint) const override;

    void FlipSurface()
    {
        double a = fA, b = fB;
        katrin::KThreeVector n1 = fN1, n2 = fN2;

        fA = b; fB = a;
        fN1 = n2; fN2 = n1;
        Update();
    }

    void Update()
    {
        fN3 = fN1.Cross(fN2).Unit();
        fP1 = fP0 + fN1 * fA;
        fP2 = fP0 + fN2 * fB;
    }

    void SetA(double d)
    {
        fA = d;
        Update();
    }
    void SetB(double d)
    {
        fB = d;
        Update();
    }
    void SetP0(const katrin::KThreeVector& p)
    {
        fP0 = p;
        Update();
    }
    void SetN1(const katrin::KThreeVector& d)
    {
        fN1 = d.Unit();
        Update();
    }
    void SetN2(const katrin::KThreeVector& d)
    {
        fN2 = d.Unit();
        Update();
    }

    double GetA() const
    {
        return fA;
    }
    double GetB() const
    {
        return fB;
    }
    const katrin::KThreeVector& GetP0() const
    {
        return fP0;
    }
    const katrin::KThreeVector& GetN1() const
    {
        return fN1;
    }
    const katrin::KThreeVector& GetN2() const
    {
        return fN2;
    }
    const katrin::KThreeVector GetN3() const
    {
        return fN3;
    }
    const katrin::KThreeVector GetP1() const
    {
        return fP1;
    }
    const katrin::KThreeVector GetP2() const
    {
        return fP2;
    }

    virtual bool ContainsPoint(const katrin::KThreeVector& aPoint) const;
    double DistanceTo(const katrin::KThreeVector& aPoint, katrin::KThreeVector& nearestPoint);

  protected:
    static bool SameSide(const katrin::KThreeVector& point, const katrin::KThreeVector& A,
                       const katrin::KThreeVector& B, const katrin::KThreeVector& C);

  protected:
    double fA, fB;
    katrin::KThreeVector fP0, fP1, fP2;
    katrin::KThreeVector fN1, fN2, fN3;
};
}  // namespace KGeoBag

#endif
