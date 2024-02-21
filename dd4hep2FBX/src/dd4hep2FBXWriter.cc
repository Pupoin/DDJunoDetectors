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
#include "XML/XML.h"

#include "Math/Vector3D.h"
#include "Math/Transform3D.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
typedef std::map<std::string, Handle<NamedObject>> HandleMap;
typedef std::map<std::string, DetElement> Children;

#include <iomanip>

// dd4hep2FBXWriter(Detector lcdd, bool usePrototypes)
// {
//   this->m_lcdd = lcdd;
//   this->m_UsePrototypes = usePrototypes;
// }

dd4hep2FBXWriter::dd4hep2FBXWriter(HandleMap det_map, bool usePrototypes)
{
  this->m_det_map = det_map;
  this->m_UsePrototypes = usePrototypes;
}

bool dd4hep2FBXWriter::doit(std::string outputFilename)
{
  // Check that the top-most (world) detector element has been created already
  // if (m_det.world() == nullptr)
  // {
  //   return false;
  // }
  // Detector &lcwdd = Detector::getInstance();
  // Tube waterPool(0, 10 * mm, 10 * mm, 20 * mm, 30 * mm);
  // lcwdd.fromCompact("/home/wln/DD4hep_source/DDDetectors/compact/SiD.xml");
  // DetElement ele=lcwdd.detector("WaterPool");
  // std::cout << "@@@@@ this is test: " << ele.id() << "  " << ele.name() << std::endl;
  /// Accessor to the map of sub-detectors
  // const HandleMap &m_det_map = m_lcdd.detectors();

  // std::vector<std::string> m_DetName;   //= new std::vector<std::string>(pvStore->size(), "");
  // std::vector<std::string> m_VolName;   //= new std::vector<std::string>(lvStore->size(), "");
  // std::vector<std::string> m_SolidName; // = new std::vector<std::string>(solidStore->size(), "");

  // Assign legal and unique names to each used physical volume, logical volume and solid
  for (const auto &[name, detHandle] : m_det_map)
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
  for (const auto &[name, detHandle] : m_det_map)
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


  // // Open the output file
  // if (outputFilename.length() > 0)
  // {
  //   m_File.open(outputFilename, std::ios_base::trunc);
  // }
  // else
  // {
  //   m_File.open("geometry.fbx", std::ios_base::trunc);
  // }
  // if (m_File.fail())
  // {
  //   return false;
  // }

  // // Write the FBX preamble and headers
  // writePreamble(modelCount, materialCount, geometryCount);

  return true;
}

std::vector<std::string> dd4hep2FBXWriter::assignName(std::vector<std::string> names, string originalName)
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

void dd4hep2FBXWriter::writePreamble(int modelCount, int materialCount, int geometryCount)
{
  std::time_t t = std::time(NULL);
  struct tm *now = std::localtime(&t);
  m_File << "; FBX 7.3.0 project file" << std::endl
         << "; Copyright (C) 1997-2010 Autodesk Inc. and/or its licensors." << std::endl
         << "; All rights reserved." << std::endl
         << std::endl
         << "FBXHeaderExtension:  {" << std::endl
         << "\tFBXHeaderVersion: 1003" << std::endl
         << "\tFBXVersion: 7300" << std::endl
         << "\tCreationTime: \"" << std::put_time(now, "%F %T") << ":000\"" << std::endl
         <<
      //"\tCreationTimeStamp:  {" << std::endl <<
      //"\t\tVersion: 1000" << std::endl <<
      //"\t\tYear: " << now->tm_year + 1900 << std::endl <<
      //"\t\tMonth: " << now->tm_mon + 1 << std::endl <<
      //"\t\tDay: " << now->tm_mday << std::endl <<
      //"\t\tHour: " << now->tm_hour << std::endl <<
      //"\t\tMinute: " << now->tm_min << std::endl <<
      //"\t\tSecond: " << now->tm_sec << std::endl <<
      //"\t\tMillisecond: 0" << std::endl <<
      //"\t}" << std::endl <<
      "\tCreator: \"FBX SDK/FBX Plugins version 2013.3\"" << std::endl
         << "\tSceneInfo: \"SceneInfo::GlobalInfo\", \"UserData\" {" << std::endl
         << "\t\tType: \"UserData\"" << std::endl
         << "\t\tVersion: 100" << std::endl
         << "\t\tMetaData:  {" << std::endl
         << "\t\t\tVersion: 100" << std::endl
         << "\t\t\tTitle: \"Belle II Detector\"" << std::endl
         << "\t\t\tSubject: \"Detector Geometry Model\"" << std::endl
         << "\t\t\tAuthor: \"Belle II Collaboration\"" << std::endl
         << "\t\t\tKeywords: \"\"" << std::endl
         << "\t\t\tRevision: \"\"" << std::endl
         << "\t\t\tComment: \"\"" << std::endl
         << "\t\t}" << std::endl
         << "\t}" << std::endl
         << "}" << std::endl
         << std::endl;

  m_File << "GlobalSettings:  {" << std::endl
         << "\tVersion: 1000" << std::endl
         << "\tProperties70:  {" << std::endl
         << "\t\tP: \"UpAxis\", \"int\", \"Integer\", \"\",1" << std::endl
         << "\t\tP: \"UpAxisSign\", \"int\", \"Integer\", \"\",1" << std::endl
         << "\t\tP: \"FrontAxis\", \"int\", \"Integer\", \"\",2" << std::endl
         << "\t\tP: \"FrontAxisSign\", \"int\", \"Integer\", \"\",1" << std::endl
         << "\t\tP: \"CoordAxis\", \"int\", \"Integer\", \"\",0" << std::endl
         << "\t\tP: \"CoordAxisSign\", \"int\", \"Integer\", \"\",1" << std::endl
         << "\t\tP: \"OriginalUpAxis\", \"int\", \"Integer\", \"\",1" << std::endl
         << "\t\tP: \"OriginalUpAxisSign\", \"int\", \"Integer\", \"\",1" << std::endl
         << "\t\tP: \"UnitScaleFactor\", \"double\", \"Number\", \"\",1" << std::endl
         << "\t\tP: \"OriginalUnitScaleFactor\", \"double\", \"Number\", \"\",1" << std::endl
         << "\t\tP: \"AmbientColor\", \"ColorRGB\", \"Color\", \"\",1,1,1" << std::endl
         << "\t\tP: \"DefaultCamera\", \"KString\", \"\", \"\", \"Producer Perspective\"" << std::endl
         << "\t\tP: \"TimeMode\", \"enum\", \"\", \"\",0" << std::endl
         << "\t\tP: \"TimeSpanStart\", \"KTime\", \"Time\", \"\",0" << std::endl
         << "\t\tP: \"TimeSpanStop\", \"KTime\", \"Time\", \"\",10" << std::endl
         << "\t\tP: \"CustomFrameRate\", \"double\", \"Number\", \"\",-1" << std::endl
         << "\t}" << std::endl
         << "}" << std::endl
         << std::endl;

  m_File << "Documents:  {" << std::endl
         << "\tCount: 1" << std::endl
         << "\tDocument: 4000000000, \"\", \"Scene\" {" << std::endl
         << "\t\tProperties70:  {" << std::endl
         << "\t\t\tP: \"SourceObject\", \"object\", \"\", \"\"" << std::endl
         << "\t\t\tP: \"ActiveAnimStackName\", \"KString\", \"\", \"\", \"\"" << std::endl
         << "\t\t}" << std::endl
         << "\t\tRootNode: 0" << std::endl
         << "\t}" << std::endl
         << "}" << std::endl
         << std::endl;

  m_File << "References:  {" << std::endl
         << "}" << std::endl
         << std::endl;

  m_File << "Definitions:  {" << std::endl
         << "\tVersion: 100" << std::endl
         << "\tCount: 4" << std::endl
         << "\tObjectType: \"GlobalSettings\" {" << std::endl
         << "\t\tCount: 1" << std::endl
         << "\t}" << std::endl;
  m_File << "\tObjectType: \"Model\" {" << std::endl
         << "\t\tCount: " << modelCount << std::endl
         << "\t\tPropertyTemplate: \"FbxNode\" {" << std::endl
         << "\t\t\tProperties70:  {" << std::endl
         << "\t\t\t\tP: \"QuaternionInterpolate\", \"enum\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationOffset\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"RotationPivot\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"ScalingOffset\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"ScalingPivot\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"TranslationActive\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"TranslationMin\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"TranslationMax\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"TranslationMinX\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"TranslationMinY\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"TranslationMinZ\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"TranslationMaxX\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"TranslationMaxY\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"TranslationMaxZ\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationOrder\", \"enum\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationSpaceForLimitOnly\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationStiffnessX\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationStiffnessY\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationStiffnessZ\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"AxisLen\", \"double\", \"Number\", \"\",10" << std::endl
         << "\t\t\t\tP: \"PreRotation\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"PostRotation\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"RotationActive\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationMin\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"RotationMax\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"RotationMinX\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationMinY\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationMinZ\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationMaxX\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationMaxY\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"RotationMaxZ\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"InheritType\", \"enum\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"ScalingActive\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"ScalingMin\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"ScalingMax\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"ScalingMinX\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"ScalingMinY\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"ScalingMinZ\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"ScalingMaxX\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"ScalingMaxY\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"ScalingMaxZ\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"GeometricTranslation\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"GeometricRotation\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"GeometricScaling\", \"Vector3D\", \"Vector\", \"\",1,1,1" << std::endl
         << "\t\t\t\tP: \"MinDampRangeX\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MinDampRangeY\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MinDampRangeZ\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MaxDampRangeX\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MaxDampRangeY\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MaxDampRangeZ\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MinDampStrengthX\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MinDampStrengthY\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MinDampStrengthZ\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MaxDampStrengthX\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MaxDampStrengthY\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"MaxDampStrengthZ\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"PreferedAngleX\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"PreferedAngleY\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"PreferedAngleZ\", \"double\", \"Number\", \"\",0" << std::endl
         << "\t\t\t\tP: \"LookAtProperty\", \"object\", \"\", \"\"" << std::endl
         << "\t\t\t\tP: \"UpVectorProperty\", \"object\", \"\", \"\"" << std::endl
         << "\t\t\t\tP: \"Show\", \"bool\", \"\", \"\",1" << std::endl
         << "\t\t\t\tP: \"NegativePercentShapeSupport\", \"bool\", \"\", \"\",1" << std::endl
         << "\t\t\t\tP: \"DefaultAttributeIndex\", \"int\", \"Integer\", \"\",0" << std::endl
         << "\t\t\t\tP: \"Freeze\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"LODBox\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"Lcl Translation\", \"Lcl Translation\", \"\", \"A\",0,0,0" << std::endl
         << "\t\t\t\tP: \"Lcl Rotation\", \"Lcl Rotation\", \"\", \"A\",0,0,0" << std::endl
         << "\t\t\t\tP: \"Lcl Scaling\", \"Lcl Scaling\", \"\", \"A\",1,1,1" << std::endl
         << "\t\t\t\tP: \"Visibility\", \"Visibility\", \"\", \"A\",1" << std::endl
         << "\t\t\t\tP: \"Visibility Inheritance\", \"Visibility Inheritance\", \"\", \"\",1" << std::endl
         << "\t\t\t}" << std::endl
         << "\t\t}" << std::endl
         << "\t}" << std::endl;
  m_File << "\tObjectType: \"Material\" {" << std::endl
         << "\t\tCount: " << materialCount << std::endl
         << "\t\tPropertyTemplate: \"FbxSurfacePhong\" {" << std::endl
         << "\t\t\tProperties70:  {" << std::endl
         << "\t\t\t\tP: \"ShadingModel\", \"KString\", \"\", \"\", \"Phong\"" << std::endl
         << "\t\t\t\tP: \"MultiLayer\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"EmissiveColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl
         << "\t\t\t\tP: \"EmissiveFactor\", \"double\", \"Number\", \"A\",0" << std::endl
         << "\t\t\t\tP: \"AmbientColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl
         << "\t\t\t\tP: \"AmbientFactor\", \"double\", \"Number\", \"A\",0" << std::endl
         << "\t\t\t\tP: \"DiffuseColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl
         << "\t\t\t\tP: \"DiffuseFactor\", \"double\", \"Number\", \"A\",1" << std::endl
         << "\t\t\t\tP: \"Bump\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"NormalMap\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"BumpFactor\", \"double\", \"Number\", \"\",1" << std::endl
         << "\t\t\t\tP: \"TransparentColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl
         << "\t\t\t\tP: \"TransparencyFactor\", \"double\", \"Number\", \"A\",0" << std::endl
         << "\t\t\t\tP: \"DisplacementColor\", \"ColorRGB\", \"Color\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"DisplacementFactor\", \"double\", \"Number\", \"\",1" << std::endl
         << "\t\t\t\tP: \"VectorDisplacementColor\", \"ColorRGB\", \"Color\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"VectorDisplacementFactor\", \"double\", \"Number\", \"\",1" << std::endl
         << "\t\t\t\tP: \"SpecularColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl
         << "\t\t\t\tP: \"SpecularFactor\", \"double\", \"Number\", \"A\",0" << std::endl
         << "\t\t\t\tP: \"ShininessExponent\", \"double\", \"Number\", \"A\",20" << std::endl
         << "\t\t\t\tP: \"ReflectionColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl
         << "\t\t\t\tP: \"ReflectionFactor\", \"double\", \"Number\", \"A\",0" << std::endl
         << "\t\t\t}" << std::endl
         << "\t\t}" << std::endl
         << "\t}" << std::endl;
  /*
  m_File << "\tObjectType: \"Material\" {" << std::endl <<
            "\t\tCount: " << materialCount << std::endl <<
            "\t\tPropertyTemplate: \"FbxSurfaceLambert\" {" << std::endl <<
            "\t\t\tProperties70:  {" << std::endl <<
            "\t\t\t\tP: \"ShadingModel\", \"KString\", \"\", \"\", \"Lambet\"" << std::endl <<
            "\t\t\t\tP: \"MultiLayer\", \"bool\", \"\", \"\",0" << std::endl <<
            "\t\t\t\tP: \"EmissiveColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl <<
            "\t\t\t\tP: \"EmissiveFactor\", \"double\", \"Number\", \"A\",0" << std::endl <<
            "\t\t\t\tP: \"AmbientColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl <<
            "\t\t\t\tP: \"AmbientFactor\", \"double\", \"Number\", \"A\",0" << std::endl <<
            "\t\t\t\tP: \"DiffuseColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl <<
            "\t\t\t\tP: \"DiffuseFactor\", \"double\", \"Number\", \"A\",1" << std::endl <<
            "\t\t\t\tP: \"Bump\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl <<
            "\t\t\t\tP: \"NormalMap\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl <<
            "\t\t\t\tP: \"BumpFactor\", \"double\", \"Number\", \"\",1" << std::endl <<
            "\t\t\t\tP: \"TransparentColor\", \"ColorRGB\", \"Color\", \"A\",0,0,0" << std::endl <<
            "\t\t\t\tP: \"TransparencyFactor\", \"double\", \"Number\", \"A\",0" << std::endl <<
            "\t\t\t\tP: \"DisplacementColor\", \"ColorRGB\", \"Color\", \"\",0,0,0" << std::endl <<
            "\t\t\t\tP: \"DisplacementFactor\", \"double\", \"Number\", \"\",1" << std::endl <<
            "\t\t\t\tP: \"VectorDisplacementColor\", \"ColorRGB\", \"Color\", \"\",0,0,0" << std::endl <<
            "\t\t\t\tP: \"VectorDisplacementFactor\", \"double\", \"Number\", \"\",1" << std::endl <<
            "\t\t\t}" << std::endl <<
            "\t\t}" << std::endl <<
            "\t}" << std::endl;
  */
  m_File << "\tObjectType: \"Geometry\" {" << std::endl
         << "\t\tCount: " << geometryCount << std::endl
         << "\t\tPropertyTemplate: \"FbxMesh\" {" << std::endl
         << "\t\t\tProperties70:  {" << std::endl
         << "\t\t\t\tP: \"Color\", \"ColorRGB\", \"Color\", \"\",1,1,1" << std::endl
         << "\t\t\t\tP: \"BBoxMin\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"BBoxMax\", \"Vector3D\", \"Vector\", \"\",0,0,0" << std::endl
         << "\t\t\t\tP: \"Primary Visibility\", \"bool\", \"\", \"\",1" << std::endl
         << "\t\t\t\tP: \"Casts Shadows\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t\tP: \"Receive Shadows\", \"bool\", \"\", \"\",0" << std::endl
         << "\t\t\t}" << std::endl
         << "\t\t}" << std::endl
         << "\t}" << std::endl
         << "}" << std::endl
         << std::endl;
}
