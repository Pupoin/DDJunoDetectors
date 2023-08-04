#include <map>
#include <stdio.h>
#include "PMT.h"
// ROOT include files
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"
#include "Math/Vector3D.h"
#include "Math/Transform3D.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

// PMT::PMT(int idx, Position pos, XYZAngles rotation, PMTType type)
// {
//     this->_idx = idx;
//     this->_pmtPosition = pos;
//     this->_pmtRotation = rotation;
//     this->_pmtType = type;
// }

PMT::PMT(Position pos, XYZAngles rotation, PMTType type)
{
    this->_pmtPosition = pos;
    this->_pmtRotation = rotation;
    this->_pmtType = type;
}

PMT::~PMT()
{
}

void PMT::setPMTType(PMTType type)
{
    this->_pmtType = type;
}

PMTType PMT::getPMTType()
{
    return _pmtType;
}

Position PMT::getPMTPosition()
{
    return _pmtPosition;
}
XYZAngles PMT::getPMTRotation()
{
    return _pmtRotation;
}