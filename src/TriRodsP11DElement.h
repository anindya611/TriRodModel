#ifndef TRI_RODS_P11DELEMENT_H
#define TRI_RODS_P11DELEMENT_H

#include <P12DElement.h>
#include <ShapesEvaluated.h>

namespace TriRods
{
  //! Linear 1D elements
  
  template<int nFields>
    class P11DElement: public Element
    {
    public:
      // Constructor
      // \param[in] i1, i2 Global numbering for nodes
      inline P11DElement(int i1, int i2)
	:Element()
      {
	  SegGeom = new Segment<1>(i1, i2);
	  const Quadrature* qrule = Line_1::Bulk;
	  Linear<1> LinShp;
	  ShapesEvaluated SE(qrule, &LinShp, SegGeom);
	  AddBasisFunctions(SE);
	  for(int i=0; i<nFields; ++i)
	    AppendField(i);
	}

      // Copy constructor
      inline P11DElement(const P11DElement<nFields>& Obj)
	:Element(Obj)
      { SegGeom = Obj.SegGeom->Clone(); }

      // Destructor
      inline ~P11DElement()
      { delete SegGeom; }

      // Cloning
      inline virtual P11DElement<nFields>* Clone() const override
      { return new P11DElement(*this); }

      // Returns the element geometry
      const Segment<1>& GetElementGeometry() const
      { return *SegGeom; }

    protected:
      int GetFieldIndex(int fields) const override
      { return 0; }

    private:
      Segment<1>* SegGeom;
    };

  // Local to global map is identical to P12DElement
  typedef StandardP12DMap StandardP11DMap;
}

#endif
