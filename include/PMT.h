#include <map>
#include <stdio.h>
// ROOT include files
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "Math/Vector3D.h"
#include "Math/Transform3D.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;


enum class PMTType {NotSet=0, HZC, Hamamatsu, NNVT, HighQENNVT};

class PMT
{
  public:


    PMT()=default;
    // PMT(int idx, Position pos, XYZAngles pmtRotation, PMTType pmttype=PMTType::NotSet);

    
    PMT(Position pos, XYZAngles pmtRotation, PMTType pmttype=PMTType::NotSet);

    PMTType getPMTType();
    void setPMTType(PMTType type);

    Position getPMTPosition();
    XYZAngles getPMTRotation();


    virtual ~PMT();




  private:
    Position  _pmtPosition=Position(0,0,0);
    XYZAngles _pmtRotation=XYZAngles(0,0,0);
    PMTType   _pmtType=PMTType::NotSet;
    // int _idx=0;



};


