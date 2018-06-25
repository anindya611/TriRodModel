#ifndef TRI_RODS_OPS_H
#define TRI_RODS_OPS_H

#include <TriRodsConfigurationNew1.h>
#include <TriRodsProperties.h>
#include <TriRodsOpsWorkspace.h>
#include <ElementalOperation.h>

namespace TriRods
{
  class Ops: public DResidue
  {
  public:
    // Constructor
    inline Ops(const Element* elm,
	       const std::vector<int> fieldnums,
	       ConfigAccessStruct str,
	       const Material& mat,
	       OpsWorkspace* wrkspc)
      :DResidue(), Elm(elm), fields(fieldnums), ConfigStr(str),
      Mat(&mat), WrkSpc(wrkspc)
    { assert(static_cast<int>(fields.size())==11 &&
	     "TriRods::Ops::Ops has expected 11 fields"); }

    // Destructor, does nothing
    inline virtual ~Ops() {}

    // Copy constructor
    inline Ops(const Ops& Obj)
      :DResidue(Obj),
      Elm(Obj.Elm), fields(Obj.fields),
      ConfigStr(Obj.ConfigStr), Mat(Obj.Mat),
      WrkSpc(Obj.WrkSpc) {}

    // Cloning
    inline virtual Ops* Clone() const override
    { return new Ops(*this); }

    // Disable assignment
    Ops& operator=(const Ops&) = delete;

    // Returns the element being used
    inline const Element* GetElement() const
    { return Elm; }

    // Returns the material being used
    inline const Material* GetMaterial() const
    { return Mat; }

    // Return the fields used
    inline const std::vector<int>& GetField() const override
    { return fields; }

    // Returns the number of dofs for a given field
    inline int GetFieldDof(int fieldnumber) const override
    { return Elm->GetDof(fields[fieldnumber]); }

    // Computes the energy functional
    double GetEnergy(const void* Config) const;

    // Residual calculation
    void GetVal(const void* Config,
		std::vector<std::vector<double>>* resval) const;

    // Residual and stiffness calculation
    void GetDVal(const void* Config,
		 std::vector<std::vector<double>>* resval,
		 std::vector<std::vector<std::vector<std::vector<double>>>>* dresval) const;

    // Consistency test for this operation at the given calculation
    bool ConsistencyTest(const void* Cpnfig,
			 const double pertEPS, const double tolEPS) const;
    
  private:
    const Element* Elm; // Element for this operation
    std::vector<int> fields; // Field numbers for this operation
    const ConfigAccessStruct ConfigStr; // Details to access configuration for this
    const Material* Mat; // Material defining the rod/shell
    OpsWorkspace* WrkSpc; // Workspace to be used for calculations
  };
}

#endif
