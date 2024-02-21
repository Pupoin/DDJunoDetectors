/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2016 - Belle II Collaboration                             *
 *                                                                        *
 * Author: ZhaoYang Yuan                                                   *
 * 2024/02/20                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

#include "dd4hep2FBXWriter.h"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Detector.h"
#include <DD4hep/Handle.h>
#include "XML/Layering.h"
#include "DD4hep/Volumes.h"
#include "XML/XML.h"

#include "Math/Vector3D.h"
#include "Math/Transform3D.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

typedef std::map<std::string, Handle<NamedObject>> HandleMap;
typedef std::map<std::string, DetElement> Children;

#include <iomanip>


int main()
{
  Detector &lcdd = Detector::getInstance();
  // Tube waterPool(0, 10 * mm, 10 * mm, 20 * mm, 30 * mm);
  lcdd.fromCompact("/home/wln/DD4hep_source/DDDetectors/compact/SiD.xml");
  const HandleMap det_map = lcdd.detectors();

  dd4hep2FBXWriter test(det_map, true);
  // dd4hep2FBXWriter test(lcdd, true);
  test.doit("test.fbx");



  return 0 ;
}

/*
std::vector<std::string> assignName(std::vector<std::string> names, string originalName)
{
  // Replace problematic characters with underscore
  for (char c : " .,:;?'\"*+-=|^!/@#$\\%{}[]()<>")
    std::replace(originalName.begin(), originalName.end(), c, '_');
  //
  for (size_t i = 0; i < names.size(); i++)
  {
    if (originalName == names[i])
    {
      names[i] = originalName + "_" + to_string(i);
    }
  }
  return names;
}

int main()
{
  Detector &lcwdd = Detector::getInstance();
  Tube waterPool(0, 10 * mm, 10 * mm, 20 * mm, 30 * mm);
  lcwdd.fromCompact("/home/wln/DD4hep_source/DDDetectors/compact/SiD.xml");
  // DetElement ele=lcwdd.detector("WaterPool");
  // std::cout << "@@@@@ this is test: " << ele.id() << "  " << ele.name() << std::endl;
  /// Accessor to the map of sub-detectors
  const HandleMap det_map = lcwdd.detectors();

  std::vector<std::string> m_DetName;   //= new std::vector<std::string>(pvStore->size(), "");
  std::vector<std::string> m_VolName;   //= new std::vector<std::string>(lvStore->size(), "");
  std::vector<std::string> m_SolidName; // = new std::vector<std::string>(solidStore->size(), "");

  // Assign legal and unique names to each used physical volume, logical volume and solid
  for (const auto &[name, detHandle] : det_map)
  {
    // get all of the names
    DetElement subdet = DetElement(detHandle);
    std::cout << __LINE__ << "name: " << name << ", detHandle: " << subdet.name() << " id " << subdet.id() << "\n";
    Volume subdetVol = subdet.volume();
    Solid subdetSolid = subdet.solid();
    std::cout << __LINE__ << "name: " << subdetVol.name() << std::endl;   //<< " id : " << subdetVol.i<< " id " << subdet.id() << "\n";
    std::cout << __LINE__ << "name: " << subdetSolid.name() << std::endl; //<< " id : " << subdetVol.i<< " id " << subdet.id() << "\n";

    // get the child except itselef( top-level here)
    // const Children&  detchd = subdet.children();
    // for (const auto &[chdname, chd] : detchd){
    //   DetElement subdetofdet = DetElement(chd);
    //   std::cout << __LINE__<<"name: " << name << ", detHandle: " << subdetofdet.name() << " id " << subdetofdet.id() << "\n";
    // }

    m_DetName.push_back(subdet.name());
    m_VolName.push_back(subdetVol.name());
    m_SolidName.push_back(subdetSolid.name());
  }
  // Assign new name if duplicate
  // Compact          ERROR ++ FAILED    to convert subdetector:
  // PMT_type1: dd4hep: Attempt to add an already existing object:PMT_type1.
  for (const auto &[name, detHandle] : det_map)
  {
    DetElement subdet = DetElement(detHandle);
    Volume subdetVol = subdet.volume();
    Solid subdetSolid = subdet.solid();

    string detName = subdet.name();
    string detVolName = subdetVol.name();
    string detSolidName = subdetSolid.name();
    
    m_DetName= assignName(m_DetName, detName);
    m_VolName= assignName(m_VolName, detName);
    m_SolidName= assignName(m_SolidName, detName);

  }

  // Count the number of references to each physical volume and logical volume and solid
  // so that these values can be placed in the FBX file's Definitions{} section.
  unsigned int geometryCount = m_DetName.size();
  unsigned int materialCount = m_VolName.size();
  unsigned int modelCount = m_SolidName.size();

  

  return 0;
}
*/