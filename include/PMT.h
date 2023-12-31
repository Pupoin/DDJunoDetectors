#include <map>
#include <stdio.h>
// ROOT include files
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "Math/Vector3D.h"
#include "Math/Transform3D.h"
#include "Math/AxisAngle.h"
#include "Math/RotationZ.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
using namespace ROOT::Math;


enum class PMTType {NotSet=0, Hamamatsu, NNVT, HighQENNVT, HZC};

class PMT
{
  public:


    PMT()=default;
    // PMT(int idx, Position pos, XYZAngles pmtRotation, PMTType pmttype=PMTType::NotSet);

    
    PMT(Position pos, AxisAngle pmtRotation, PMTType pmttype=PMTType::NotSet);

    PMTType getPMTType();
    void setPMTType(PMTType type);

    Position getPMTPosition();
    AxisAngle getPMTRotation();


    virtual ~PMT();




  private:
    Position  _pmtPosition=Position(0,0,0);
    AxisAngle _pmtRotation;
    PMTType   _pmtType=PMTType::NotSet;
    // int _idx=0;



};


