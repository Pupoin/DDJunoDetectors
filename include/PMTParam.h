#include <map>
#include <vector>
#include <stdio.h>
// #include "PMT.h"
// DD4hep include files
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



class PMTParam 
{
  public:
    PMTParam()=default;


    PMTParam(std::string positionFilePath, PMTType defaultPMTType);


    /// @param positionFilePath PMT position file path
    /// @param typeFilePath  PMT type file path
    PMTParam(std::string positionFilePath, std::string typeFilePath);


    map<int,PMT> getPMTParam();
    virtual ~PMTParam();

  private:
    void readFileGetPosition(int paramNumPerLine);
    void readFileGetType(int paramNumPerLine);


    
    private:

      std::string _positionFilePath="", _typeFilePath="";
      PMTType _defaultPMTType=PMTType::NotSet;
      map<int,PMT> _pmtParam{};
};

