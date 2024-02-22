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
#include "HepPolyhedron.h"
#include "G4Polyhedron.hh"

#include "Math/Vector3D.h"
#include "Math/Transform3D.h"
#include <iomanip>

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
typedef std::map<std::string, Handle<NamedObject>> HandleMap;
typedef std::map<std::string, DetElement> Children;

dd4hep2FBXWriter::dd4hep2FBXWriter(string filePath, bool usePrototypes)
{
  this->m_filePath = filePath;
  this->m_UsePrototypes = usePrototypes;
}

void dd4hep2FBXWriter::getAllChildren(DetElement det)
{
  m_childrenDet.push_back(det);
  m_childrenVol.push_back(det.volume());
  m_childrenSolid.push_back(det.solid());

  m_DetName.push_back(det.name());
  m_VolName.push_back(det.volume().name());
  m_SolidName.push_back(det.solid().name());

  // std::cout << " is assembly " << det.volume().isAssembly() << std::endl;
  // std::cout << " in getname() " << det.name() << std::endl;
  if (det.children().size() != 0)
  {
    for (const auto &[name, detchild] : det.children())
    {
      getAllChildren(detchild);
    }
  }
}

bool dd4hep2FBXWriter::doit(std::string outputFilename)
{

  Detector &m_lcdd = Detector::getInstance();
  // Tube waterPool(0, 10 * mm, 10 * mm, 20 * mm, 30 * mm);
  m_lcdd.fromCompact(m_filePath);
  // DetElement ele=lcwdd.detector("WaterPool");
  // std::cout << "@@@@@ this is test: " << ele.id() << "  " << ele.name() << std::endl;
  /// Accessor to the map of sub-detectors
  // m_det_map = m_lcdd.detectors();
  m_world = m_lcdd.world();
  // Check that the top-most (world) detector element has been created already
  if (!m_lcdd.world())
  {
    return false;
  }

  // Assign legal and unique names to each used physical volume, logical volume and solid
  getAllChildren(m_world);
  // Assign new name if duplicate
  // Compact          ERROR ++ FAILED    to convert subdetector:
  // PMT_type1: dd4hep: Attempt to add an already existing object:PMT_type1.
  // for (const auto subdet : m_childrenDet)
  for (size_t i = 0; i < m_childrenDet.size(); i++)
  {
    DetElement subdet = m_childrenDet[i];
    // DetElement subdet = DetElement(mdet);
    // std::cout << __LINE__ << " name:" << subdet.name() << " children:" << subdet.children().size()
    //           << " isassembly: " << subdet.volume().isAssembly() << std::endl;

    Volume subdetVol = subdet.volume();
    Solid subdetSolid = subdet.solid();

    // auto c=subdetVol.IsReplicated();

    string detName = subdet.name();
    string detVolName = subdetVol.name();
    string detSolidName = subdetSolid.name();

    m_DetName = assignName(m_DetName, detName, i);
    m_VolName = assignName(m_VolName, detName, i);
    m_SolidName = assignName(m_SolidName, detName, i);
  }

  // Count the number of references to each physical volume and logical volume and solid
  // so that these values can be placed in the FBX file's Definitions{} section.
  // countEntities(m_world);
  unsigned int geometryCount = m_DetName.size();
  unsigned int materialCount = m_VolName.size();
  unsigned int modelCount = m_SolidName.size();
  // Open the output file
  if (outputFilename.length() > 0)
  {
    m_File.open(outputFilename, std::ios_base::trunc);
  }
  else
  {
    m_File.open("geometry.fbx", std::ios_base::trunc);
  }
  if (m_File.fail())
  {
    return false;
  }

  // Write the FBX preamble and headers
  writePreamble(modelCount, materialCount, geometryCount);

  // Write all solids as Geometry nodes (replicas are written later).
  // Write all logical volumes as Material nodes (color information).
  // Write all physical and logical volumes as Model nodes (with replica-solids treated here).
  // m_PVID = new std::vector<unsigned long long>(pvStore->size(), 0x0000010000000000LL);
  // m_LVID = new std::vector<unsigned long long>(lvStore->size(), 0x000000C000000000LL);
  m_SolidID = new std::vector<unsigned long long>(m_childrenSolid.size(), 0x0000008000000000LL);
  // m_MatID = new std::vector<unsigned long long>(lvStore->size(), 0x0000004000000000LL);
  // m_Visible = new std::vector<bool>(lvStore->size(), false);

  m_File << "Objects:  {" << std::endl;
  for (unsigned int solidIndex = 0; solidIndex <m_childrenSolid.size(); ++solidIndex)
  {
    (*m_SolidID)[solidIndex] += 0x0000000001000000LL * solidIndex;
    if (m_SolidName[solidIndex].length() > 0)
    {
      // for (unsigned int solidCount = 0; solidCount <= m_childrenSolid[solidIndex]; ++solidCount) { // note lower and upper limits!
      std::cout << __LINE__ << " " << m_childrenSolid.size() << " " << m_SolidID->size()
                << m_SolidName[solidIndex] << " "
                << " index " << solidIndex << " "
                << (*m_SolidID)[solidIndex] << std::endl;
      // cout << m_childrenSolid[solidIndex].to_string() << " " << 
      writeGeometryNode(m_childrenSolid[solidIndex], m_SolidName[solidIndex], (*m_SolidID)[solidIndex]);
      // }
      break;
    }
  }

  cout << __LINE__ << " " << m_DetName.size() << std::endl;

  return true;
}
/*
void dd4hep2FBXWriter::countEntities(DetElement world)
{
  G4VPhysicalVolume* physVol;
  // Descend to the leaves of the tree
  G4LogicalVolume* logVol = physVol->GetLogicalVolume();
  for (int daughter = 0; daughter < logVol->GetNoDaughters(); ++daughter) {
    G4VPhysicalVolume* physVolDaughter = logVol->GetDaughter(daughter);
    for (int j = 0; j < physVolDaughter->GetMultiplicity(); ++j) {
      countEntities(physVolDaughter);
    }
  }
  // Count replicas and duplicates of each physical and logical volume as well as the unique
  // versions of replicated solids as we ascend the recursive tree
  G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  G4SolidStore* solidStore = G4SolidStore::GetInstance();
  G4VSolid* solid = logVol->GetSolid();
  int pvIndex = std::find(pvStore->begin(), pvStore->end(), physVol) - pvStore->begin();
  int lvIndex = std::find(lvStore->begin(), lvStore->end(), logVol) - lvStore->begin();
  int solidIndex = std::find(solidStore->begin(), solidStore->end(), solid) - solidStore->begin();
  if (physVol->IsReplicated())
  {
    EAxis axis;
    G4int nReplicas;
    G4double width;
    G4double offset;
    G4bool consuming;
    physVol->GetReplicationData(axis, nReplicas, width, offset, consuming);
    G4VPVParameterisation *physParameterisation = physVol->GetParameterisation();
    if (physParameterisation)
    { // parameterised volume
      G4VSolid *solidReplica = physParameterisation->ComputeSolid(0, physVol);
      physParameterisation->ComputeTransformation(0, physVol);
      solidReplica->ComputeDimensions(physParameterisation, 0, physVol);
      if (!(*solidReplica == *solid))
        (*m_SolidReplicas)[solidIndex]++;
      if (m_UsePrototypes && (*solidReplica == *solid))
      {
        if ((*m_LVReplicas)[lvIndex] > 0)
          (*m_LVUnique)[lvIndex] = false;
        (*m_LVReplicas)[lvIndex] = 1;
      }
      else
      {
        (*m_LVReplicas)[lvIndex]++;
      }
      (*m_PVReplicas)[pvIndex]++;
    }
    else
    { // plain replicated volume
      if ((axis == kRho) && (solid->GetEntityType() == "G4Tubs"))
        (*m_SolidReplicas)[solidIndex]++;
      if (m_UsePrototypes && !((axis == kRho) && (solid->GetEntityType() == "G4Tubs")))
      {
        (*m_LVReplicas)[lvIndex] = 1;
      }
      else
      {
        (*m_LVReplicas)[lvIndex]++;
      }
      (*m_PVReplicas)[pvIndex]++;
    }
  }
  else
  {
    if ((*m_LVCount)[lvIndex] > 0)
      (*m_LVUnique)[lvIndex] = false;
    if (m_UsePrototypes)
    {
      (*m_PVCount)[pvIndex] = 1;
      (*m_LVCount)[lvIndex] = 1;
    }
    else
    {
      (*m_PVCount)[pvIndex]++;
      (*m_LVCount)[lvIndex]++;
    }
  }
}
*/

void getPolyhedron(Solid solid)
{
    // bool solidtype=solid->IsCylType();
  string solidtype=solid.type();
  if (solidtype == "TGeoBBox")
  {
    auto *ddda = new HepPolyhedronBox(Box(solid).x(), Box(solid).y(), Box(solid).z());
    /* code */

  }

  // return ;
}

void dd4hep2FBXWriter::writeGeometryNode(Solid solid, const std::string solidName, unsigned long long solidID)
{
  std::cout <<__LINE__ << " solid.type:" << solid.type() << /*" IsComposite:" << solid.IsComposite() << */ std::endl;
  string solidtype=solid.type();
  if ((solidtype == "TGeoIntersection") ||
      (solidtype == "TGeoUnion") ||
      (solidtype == "TGeoSubtraction") ||
      (solidtype == "TGeoBoolNode")) {
    // HepPolyhedron* polyhedron = getBooleanSolidPolyhedron(solid);
    // G4Polyhedron* g4polyhedron = new G4Polyhedron(*polyhedron);
    // writePolyhedron(solid, g4polyhedron, solidName, solidID);
    // delete polyhedron;
    // delete g4polyhedron;
  } else {
    auto a=Polyhedra(solid);
    std:: cout << __LINE__ << " name:" << solid.name() << " type:" << solid.type()<< std::endl;
    std:: cout << __LINE__
      << "solid.tostring: " << solid.toString()
      << " solid->x()" << Box(solid).x() 
      // << " solid->y()" << solid.y() 
      // << " solid->z()" << solid.z() 
      //  << " name:" << a.name() 
      // << " type:" << a.type()
      << std::endl;
    
    getPolyhedron(solid);
    // std::cout << __LINE__ << " polyhedra:" << a.numEdges() << std::endl;
    // writePolyhedron(solid, PolyhedraRegular(solid), solidName, solidID);
  }
}


/*
void dd4hep2FBXWriter::writePolyhedron(Solid solid, TGeoPolygon* polyhedron, const std::string name,
                                      unsigned long long solidID)
{
  if (polyhedron) {
    polyhedron->SetNumberOfRotationSteps(120);
    m_File << "\t; Solid " << solid->GetName() << " of type " << solid->GetEntityType() << std::endl <<
           "\tGeometry: " << solidID << ", \"Geometry::" << name << "\", \"Mesh\" {" << std::endl <<
           "\t\tVertices: *" << polyhedron->GetNoVertices() * 3 << " {" << std::endl << "\t\t\ta: ";
    std::streampos startOfLine = m_File.tellp();
    for (int j = 1; j <= polyhedron->GetNoVertices(); ++j) {
      m_File << (j == 1 ? "" : ",") <<
             polyhedron->GetVertex(j).x() << "," <<
             polyhedron->GetVertex(j).y() << "," <<
             polyhedron->GetVertex(j).z();
      if (m_File.tellp() - startOfLine > 100) {
        startOfLine = m_File.tellp();
        m_File << std::endl << "\t\t\t\t";
      }
    }
    m_File << std::endl << "\t\t}" << std::endl;

    std::vector<int> vertices;
    for (int k = 1; k <= polyhedron->GetNoFacets(); ++k) {
      G4bool notLastEdge = true;
      G4int ndx = -1, edgeFlag = 1;
      do {
        notLastEdge = polyhedron->GetNextVertexIndex(ndx, edgeFlag);
        if (notLastEdge) {
          vertices.push_back(ndx - 1);
        } else {
          vertices.push_back(-ndx);
        }
      } while (notLastEdge);
    }
    m_File << "\t\tPolygonVertexIndex: *" << vertices.size() << " {" << std::endl << "\t\t\ta: ";
    startOfLine = m_File.tellp();
    for (unsigned int j = 0; j < vertices.size(); ++j) {
      m_File << (j == 0 ? "" : ",") << vertices[j];
      if (m_File.tellp() - startOfLine > 100) {
        startOfLine = m_File.tellp();
        m_File << std::endl << "\t\t\t\t";
      }
    }
    m_File << std::endl << "\t\t}" << std::endl;

    m_File << "\t\tGeometryVersion: 124" << std::endl <<
           "\t\tLayerElementNormal: 0 {" << std::endl <<
           "\t\t\tVersion: 101" << std::endl <<
           // "\t\t\tName: \"\"" << std::endl <<
           "\t\t\tMappingInformationType: \"ByPolygonVertex\"" << std::endl <<
           "\t\t\tReferenceInformationType: \"Direct\"" << std::endl <<
           "\t\t\tNormals: *" << vertices.size() * 3 << " {" << std::endl << "\t\t\t\ta: ";
    startOfLine = m_File.tellp();
    unsigned int j = 0;
    for (int k = 1; k <= polyhedron->GetNoFacets(); ++k) {
      G4Normal3D normal = polyhedron->GetUnitNormal(k);
      do {
        m_File << (j == 0 ? "" : ",") << normal.x() << "," << normal.y() << "," << normal.z();
        if (m_File.tellp() - startOfLine > 100) {
          startOfLine = m_File.tellp();
          m_File << std::endl << "\t\t\t\t";
        }
      } while (vertices[j++] >= 0);
    }
    m_File << std::endl << "\t\t\t}" << std::endl << "\t\t}" << std::endl <<
           "\t\tLayerElementMaterial: 0 {" << std::endl <<
           "\t\t\tVersion: 101" << std::endl <<
           // "\t\t\tName: \"\"" << std::endl <<
           "\t\t\tMappingInformationType: \"AllSame\"" << std::endl <<
           "\t\t\tReferenceInformationType: \"IndexToDirect\"" << std::endl <<
           "\t\t\tMaterials: *1 {" << std::endl <<
           "\t\t\t\ta: 0" << std::endl <<
           "\t\t\t}" << std::endl <<
           "\t\t}" << std::endl <<
           "\t\tLayer: 0 {" << std::endl <<
           "\t\t\tVersion: 100" << std::endl <<
           "\t\t\tLayerElement:  {" << std::endl <<
           "\t\t\t\tType: \"LayerElementNormal\"" << std::endl <<
           "\t\t\t\tTypedIndex: 0" << std::endl <<
           "\t\t\t}" << std::endl <<
           "\t\t\tLayerElement:  {" << std::endl <<
           "\t\t\t\tType: \"LayerElementMaterial\"" << std::endl <<
           "\t\t\t\tTypedIndex: 0" << std::endl <<
           "\t\t\t}" << std::endl <<
           "\t\t}" << std::endl <<
           "\t}" << std::endl;
  } else {
    std::err << "Polyhedron representation of solid " << name << " cannot be created" < std::endl;
  }
}
*/
std::vector<std::string> dd4hep2FBXWriter::assignName(std::vector<std::string> names, string originalName, unsigned int mindex)
{
  // Replace problematic characters with underscore
  if (originalName.length() == 0)
  {
    originalName = "anonymous";
  }
  for (char c : " .,:;?'\"*+-=|^!/@#$\\%{}[]()<>")
    std::replace(originalName.begin(), originalName.end(), c, '_');
  //
  for (size_t i = 0; i < names.size(); i++)
  {
    if (i == mindex)
      continue;
    if (originalName == names[i])
    {
      names[i] = originalName + "_" + to_string(i);
      // std::cout << __LINE__ << " " << originalName << " " << names[i] << std::endl;
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
