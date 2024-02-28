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
#include "GPolyhedron.hh"
#include <TGeoBoolNode.h>
#include <TGeoScaledShape.h>
#include <TGeoCompositeShape.h>
#include "TColor.h"
#include "TGeoMatrix.h"
#include <CLHEP/Vector/Rotation.h>
#include <CLHEP/Geometry/Transform3D.h>
#include "TROOT.h"
// #include "Math/Vector3D.h"
// #include "Math/Transform3D.h"
#include <iomanip>
#include <typeinfo>
#include <CLHEP/Geometry/Normal3D.h> //#include "HepGeom::Normal3D<double>.hh"

// typedef HepGeom::Normal3D<double> G4Normal3D;

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

void dd4hep2FBXWriter::getDets(DetElement det)
{
  m_det.push_back(det);
  // m_placedVol.push_back(det.placement());
  // m_vol.push_back(det.volume());
  // m_solid.push_back(det.solid());

  // m_detName.push_back(det.name());
  // m_volName.push_back(det.volume().name());
  // m_VolChildrenNum.push_back(det.volume()->GetNdaughters());
  // det.volume()->setName("adf");
  // m_solidName.push_back(det.solid().name());
  // det.solid().setName("adf");
  // std::cout << __LINE__ << " volume: " << det.volume().name() << " " << det.volume()->GetNdaughters() << std::endl;
  // std::cout << __LINE__ << " solid: " << det.solid().name() << " " << det.solid().type() << std::endl;
  if (det.children().size() != 0)
  {
    for (const auto &[name, detchild] : det.children())
    {
      getDets(detchild);
    }
  }
}
// void dd4hep2FBXWriter::getVolSolid(TGeoNode *node) {
// vis attributes

// }
void dd4hep2FBXWriter::getVolSolid(TGeoNode *node)
{
  double nodeDaughters = node->GetNdaughters();

  TGeoVolume *aa = node->GetVolume();                 //->IsVisible() ;
  TColor *colo = gROOT->GetColor(aa->GetFillColor()); //(bb->GetFillColor());
  colo->SetAlpha(1 - (int(aa->GetTransparency()) / 100.0));

  std::cout << __LINE__ << " volname:" << node->GetName() << "  nodeDaughtersN:" << nodeDaughters
            << " volname:" << node->GetVolume()->GetName()
            << " visiable:" << aa->IsVisible()
            << " visDaughter:" << aa->IsVisDaughters()
            << " lineColor:" << aa->GetLineColor()
            << " FillColor:" << aa->GetFillColor()
            << " ,r g b " << colo->GetRed() << " " << colo->GetGreen() << " " << colo->GetBlue()
            << " alpha:" << 1 - int(aa->GetTransparency()) / 100.0
            << std::endl;

  m_placedVol.push_back(node);
  m_placedVolName.push_back(node->GetName());
  m_vol.push_back(node->GetVolume());
  m_volName.push_back(node->GetVolume()->GetName());
  m_volChildrenNum.push_back(nodeDaughters);
  m_solid.push_back(node->GetVolume()->GetShape());
  m_solidName.push_back(node->GetVolume()->GetShape()->GetName());
  m_color.push_back(colo);
  m_visible.push_back(aa->IsVisible());

  bool isVisDaughter = node->GetVolume()->IsVisDaughters();
  for (int i = 0; i < nodeDaughters; i++)
  {
    TGeoNode *daug = node->GetDaughter(i);
    if (daug->GetNdaughters() == 0)
    {
      TGeoVolume *bb = daug->GetVolume();
      TColor *color = gROOT->GetColor(bb->GetFillColor());
      color->SetAlpha(1 - (int(bb->GetTransparency()) / 100.0));

      m_placedVol.push_back(daug);
      m_placedVolName.push_back(daug->GetName());
      m_vol.push_back(daug->GetVolume());
      m_volName.push_back(daug->GetVolume()->GetName());
      m_volChildrenNum.push_back(0);
      m_solid.push_back(daug->GetVolume()->GetShape());
      m_solidName.push_back(daug->GetVolume()->GetShape()->GetName());
      m_color.push_back(color);
      if(isVisDaughter){
        m_visible.push_back(bb->IsVisible());
      }
      else{
        m_visible.push_back(false);
      }

      std::cout << __LINE__ << " vol:" << daug->GetVolume()->GetName()
                << " mother:" << daug->GetMotherVolume()->GetName()
                << " solid:" << daug->GetVolume()->GetShape()->GetName()
                << " type:" << daug->GetVolume()->GetShape()->IsA()
                << " volname:" << daug->GetVolume()->GetName()
                << " visiable:" << bb->IsVisible()
                << " visDaughter:" << bb->IsVisDaughters()
                << " lineColor:" << bb->GetLineColor()
                << " FillColor:" << bb->GetFillColor()
                << " ,r g b " << color->GetRed() << " " << color->GetGreen() << " " << color->GetBlue()
                << " alpha:" << 1 - int(bb->GetTransparency()) / 100.0
                << std::endl;
    }
    else
    {
      getVolSolid(daug);
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
  // std::cout << __LINE__ << " " <<
  // m_world.volume()->GetNode(0)->GetVolume()->GetName() << " " <<
  // m_world.volume()->GetNode(1)->GetVolume()->GetName() << std::endl;

  // Assign legal and unique names to each used physical volume, logical volume and solid
  getDets(m_world);
  getVolSolid(&*(m_world.placement()));
  std::cout << " m_placedVol size " << m_placedVol.size() << std::endl;
  for (size_t i = 0; i < m_placedVol.size(); i++)
  {
    std::cout << " $$$  phvol:" << m_placedVolName[i] << " vol:" << m_volName[i]
              << "solid:" << m_solidName[i] << " childNum:" << m_volChildrenNum[i]
              << " color:argb:" << m_color[i]->GetAlpha() << " " << m_color[i]->GetRed() << " " << m_color[i]->GetGreen() << " " << m_color[i]->GetBlue()
              << " visiable:" << m_visible[i]
              << std::endl;
  }
  // for(auto node : m_placedVol){
  //   // std::cout << node->GetVolume()->GetName() << std::endl;

  // }

  std::cout << "***********************" << std::endl;

  // Assign new name if duplicate
  for (size_t i = 0; i < m_det.size(); i++)
  {
    DetElement subdet = m_det[i];
    string detName = subdet.name();
    m_detName = assignName(m_detName, detName, i);
  }
  for (size_t i = 0; i < m_vol.size(); i++)
  {
    string detVolName = m_volName[i];
    string detSolidName = m_solidName[i];

    m_volName = assignName(m_volName, detVolName, i);
    m_solidName = assignName(m_solidName, detSolidName, i);
  }
  // Count the number of references to each physical volume and logical volume and solid
  // so that these values can be placed in the FBX file's Definitions{} section.
  // countEntities(m_world);
  unsigned int geometryCount = m_solid.size();
  unsigned int materialCount = m_vol.size();
  unsigned int modelCount = m_placedVol.size(); // wrong
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
  const unsigned long long NN = m_vol.size();
  m_PVID = new std::vector<unsigned long long>(NN, 0x0000010000000000LL);
  m_LVID = new std::vector<unsigned long long>(NN, 0x000000C000000000LL);
  m_SolidID = new std::vector<unsigned long long>(NN, 0x0000008000000000LL);
  m_MatID = new std::vector<unsigned long long>(NN, 0x0000004000000000LL);
  // m_Visible = new std::vector<bool>(m_vol.size(), false);

  m_File << "Objects:  {" << std::endl;
  // index=0 is the world,

  for (unsigned int solidIndex = 0; solidIndex < m_solid.size(); ++solidIndex)
  {
    (*m_SolidID)[solidIndex] += 0x0000000001000000LL * solidIndex;
    if (m_solidName[solidIndex].length() > 0)
    {
      std::cout << __LINE__ << " " << m_solid.size() << " " << m_SolidID->size()
                << m_solidName[solidIndex] << " "
                << " index " << solidIndex << " "
                << (*m_SolidID)[solidIndex] << std::endl;
      writeGeometryNode(m_solid[solidIndex], m_solidName[solidIndex], (*m_SolidID)[solidIndex]);
    }
  }

  // write materials
  for (unsigned int lvIndex = 0; lvIndex < m_vol.size(); ++lvIndex)
  {
    (*m_MatID)[lvIndex] += 0x0000000001000000LL * lvIndex;
    (*m_LVID)[lvIndex] += 0x0000000001000000LL * lvIndex;
    // if (m_volName[lvIndex].length() > 0) {
    // if (!(*m_LVUnique)[lvIndex])
    std::cout << __LINE__ << " materials_childVol_index: " << lvIndex << std::endl;
    writeMaterialNode(lvIndex, m_volName[lvIndex]);

    // }
  }

  for (unsigned int pvIndex = 0; pvIndex < NN; ++pvIndex)
  {
    (*m_PVID)[pvIndex] += 0x0000000001000000LL * pvIndex;
  }

  //

  // m_PVCount = new std::vector<unsigned int>(pvStore->size(), 0);
  // m_LVCount = new std::vector<unsigned int>(lvStore->size(), 0);
  // m_SolidCount = new std::vector<unsigned int>(solidStore->size(), 0);
  for (unsigned int i = 0; i < m_placedVol.size(); i++)
  {
    // std::cout << __LINE__ << " addmodels " << i << std::endl;
    addModels(m_placedVol[i], 0, i);
  }
  m_File << "}" << std::endl
         << std::endl;

  // Recursively write the connections among the solid and logical/physical volume elements
  m_PVCount = new std::vector<unsigned int>(m_detName.size(), 0);
  m_LVCount = new std::vector<unsigned int>(m_volName.size(), 0);
  m_SolidCount = new std::vector<unsigned int>(m_solidName.size(), 0);
  m_File << "Connections:  {" << std::endl;
  // addConnections(m_world, 0);
  for (unsigned int i = 0; i < m_placedVol.size(); i++)
  {
    double pvIndex = i;
    unsigned long long pvID = (*m_PVID)[pvIndex];
    std::string pvName = m_placedVolName[pvIndex];

    TGeoVolume *mother = m_placedVol[i]->GetMotherVolume();
    double index = std::find(m_vol.begin(), m_vol.end(), mother) - m_vol.begin();
    std::string mothername = m_volName[index];
    unsigned long long lvID = (*m_LVID)[index];

    writePVToParentPV(pvName, mothername, pvID, lvID);
  }
  for (unsigned int i = 0; i < m_solid.size(); i++)
  {
    double lvIndex = i;
    unsigned long long lvID = (*m_LVID)[lvIndex];
    std::string lvName = m_volName[lvIndex];
    // double pvIndex = i;
    unsigned long long solidID = (*m_PVID)[i];
    std::string solidName = m_solidName[i];

    // writeSolidToLV(lvName, solidName, (*m_Visible)[lvIndex], matID, lvID, solidID);
    // writeSolidToLV(lvName, solidName, true, matID, lvID, solidID);
  }
  for (unsigned int i = 0; i < m_vol.size(); i++)
  {
    double lvIndex = i;
    unsigned long long lvID = (*m_LVID)[lvIndex];
    std::string lvName = m_volName[lvIndex];
    double pvIndex = i;
    unsigned long long pvID = (*m_PVID)[pvIndex];
    std::string pvName = m_placedVolName[pvIndex];

    writeLVToPV(pvName, lvName, pvID, lvID);
  }

  int pvIndex = std::find(m_det.begin(), m_det.end(), m_world) - m_det.begin();
  m_File << "\t; Physical volume Model::" << m_world.name() << " to Model::RootNode" << std::endl
         << "\tC: \"OO\"," << (*m_PVID)[pvIndex] << ",0" << std::endl
         << std::endl
         << "}" << std::endl
         << std::endl;

  m_File << "Takes:  {" << std::endl
         << "\tCurrent: \"\"" << std::endl
         << "}" << std::endl;

  m_File.close();

  return true;
}

void dd4hep2FBXWriter::addConnections(DetElement physVol, int replica)
{
  // Write the PhysVolModel-parentLogVolModel connections as we descend the recursive tree.
  // If the parentLogVol is referenced at most once, use its referencing PhysVol instead.
  // G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
  // G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  // G4SolidStore* solidStore = G4SolidStore::GetInstance();
  Volume logVol = physVol.volume(); // GetLogicalVolume();
  int pvIndex = std::find(m_det.begin(), m_det.end(), physVol) - m_det.begin();
  ; //
  /// std::find(pvStore->begin(), pvStore->end(), physVol) - pvStore->begin();
  unsigned long long pvID = (*m_PVID)[pvIndex];
  unsigned int pvCount = (*m_PVCount)[pvIndex];
  std::string pvName = m_detName[pvIndex];
  int lvIndex = pvIndex; // std::find(lvStore->begin(), lvStore->end(), logVol) - lvStore->begin();
  unsigned long long lvID = (*m_LVID)[lvIndex];
  unsigned int lvCount = (*m_LVCount)[lvIndex];
  std::string lvName = m_volName[lvIndex];
  Children mdetChildren = physVol.children();
  // for (long unsigned int daughter = 0; daughter < mdetChildren.size(); ++daughter) {
  for (const auto &[name, physVolDaughter] : mdetChildren)
  {
    // DetElement physVolDaughter = mdetChildren.at(daughter).second ;// ->GetDaughter(daughter);
    int pvIndexDaughter = std::find(m_det.begin(), m_det.end(), physVolDaughter) - m_det.begin();
    unsigned long long pvIDDaughter = (*m_PVID)[pvIndexDaughter];
    unsigned int pvCountDaughter = physVolDaughter.volume()->GetNdaughters();
    // multiple phyVols in one logicalVol
    for (int j = 0; j < pvCountDaughter + 1; ++j)
    {
      //   if (m_UsePrototypes) {
      if ((replica == 0) && (j == 0) && (lvCount == 0) && (pvCountDaughter == 0))
      {
        if (mdetChildren.size() != 0)
        {
          writePVToParentPV(m_detName[pvIndexDaughter], pvName, pvIDDaughter, pvID);
        }
        else
        {
          writePVToParentLV(m_detName[pvIndexDaughter], lvName, pvIDDaughter, lvID);
        }
      }
      // } else {
      //   //writePVToParentLV(m_detName[pvIndexDaughter], lvName, pvIDDaughter+0x00010000*j+pvCountDaughter, lvID+0x00010000*replica+lvCount);
      //   writePVToParentPV(m_detName[pvIndexDaughter], pvName, pvIDDaughter + 0x00010000 * j + pvCountDaughter,
      //                     pvID + 0x00010000 * replica + pvCount);
      // }
      addConnections(physVolDaughter, j);
    }
  }

  // Write the Geometry-LogVolModel, Material-LogVolModel and PhysVolModel-LogVolModel
  // connections as we ascend the recursive tree
  Solid solid = logVol.solid();
  int solidIndex = pvIndex; // std::find(solidStore->begin(), solidStore->end(), solid) - solidStore->begin();
  unsigned long long solidID = (*m_SolidID)[solidIndex];
  unsigned long long matID = (*m_MatID)[lvIndex];
  std::string solidName = m_solidName[solidIndex];
  // if (physVol->IsReplicated()) {
  //   pvName.append("_R");
  //   pvName.append(std::to_string(replica));
  //   EAxis axis;
  //   G4int nReplicas;
  //   G4double width;
  //   G4double offset;
  //   G4bool consuming;
  //   physVol->GetReplicationData(axis, nReplicas, width, offset, consuming);
  //   physVol->SetCopyNo(replica);
  //   G4VPVParameterisation* physParameterisation = physVol->GetParameterisation();
  //   if (physParameterisation) { // parameterised volume
  //     G4VSolid* solidReplica = physParameterisation->ComputeSolid(replica, physVol);
  //     physParameterisation->ComputeTransformation(replica, physVol);
  //     solidReplica->ComputeDimensions(physParameterisation, replica, physVol);
  //     if (!(*solidReplica == *solid)) {
  //       solidName.append("_R");
  //       solidName.append(std::to_string(replica));
  //       solidID += 0x00010000 * replica;
  //     }
  //     if (m_UsePrototypes && (*solidReplica == *solid)) {
  //       if ((replica == 0) && (lvCount == 0)) {
  //         if ((*m_LVUnique)[lvIndex]) { // bypass the singleton logical volume
  //           writeSolidToPV(pvName, solidName, (*m_Visible)[lvIndex], matID, pvID, solidID);
  //         } else {
  //           writeSolidToLV(lvName, solidName, (*m_Visible)[lvIndex], matID, lvID, solidID);
  //         }
  //       }
  //     } else {
  //       lvName.append("_R");
  //       lvName.append(std::to_string(replica));
  //       if ((*m_LVUnique)[lvIndex]) { // bypass the singleton logical volume
  //         writeSolidToPV(pvName, solidName, (*m_Visible)[lvIndex], matID, pvID + 0x00010000 * replica + pvCount, solidID);
  //       } else {
  //         writeSolidToLV(lvName, solidName, (*m_Visible)[lvIndex], matID, lvID + 0x00010000 * replica + lvCount, solidID);
  //       }
  //     }
  //     if (!(*m_LVUnique)[lvIndex]) {
  //       writeLVToPV(pvName, lvName, pvID + 0x00010000 * replica + pvCount, lvID + 0x00010000 * replica + lvCount);
  //     }
  //   } else { // plain replicated volume
  //     if ((axis == kRho) && (solid->GetEntityType() == "G4Tubs")) {
  //       solidName.append("_R");
  //       solidName.append(std::to_string(replica));
  //       solidID += 0x00010000 * replica;
  //     }
  //     if (m_UsePrototypes && !((axis == kRho) && (solid->GetEntityType() == "G4Tubs"))) {
  //       if ((replica == 0) && (lvCount == 0)) {
  //         if ((*m_LVUnique)[lvIndex]) { // bypass the singleton logical volume
  //           writeSolidToPV(pvName, solidName, (*m_Visible)[lvIndex], matID, pvID, solidID);
  //         } else {
  //           writeSolidToLV(lvName, solidName, (*m_Visible)[lvIndex], matID, lvID, solidID);
  //         }
  //       }
  //     } else {
  //       lvName.append("_R");
  //       lvName.append(std::to_string(replica));
  //       if ((*m_LVUnique)[lvIndex]) { // bypass the singleton logical volume
  //         writeSolidToPV(pvName, solidName, (*m_Visible)[lvIndex], matID, pvID + 0x00010000 * replica + pvCount, solidID);
  //       } else {
  //         writeSolidToLV(lvName, solidName, (*m_Visible)[lvIndex], matID, lvID + 0x00010000 * replica + lvCount, solidID);
  //       }
  //     }
  //     if (!(*m_LVUnique)[lvIndex]) {
  //       writeLVToPV(pvName, lvName, pvID + 0x00010000 * replica + pvCount, lvID + 0x00010000 * replica + lvCount);
  //     }
  //   }
  // } else
  {
    // if (m_UsePrototypes)
    // {
    //   if (lvCount == 0)
    //   {
    //     if (physVol.children().size() == 0)
    //     { // bypass the singleton logical volume
    //       writeSolidToPV(pvName, solidName, (*m_Visible)[lvIndex], matID, pvID, solidID);
    //     }
    //     else
    //     {
    writeSolidToLV(lvName, solidName, m_visible[lvIndex], matID, lvID, solidID);
    //     }
    //   }
    //   if (pvCount == 0)
    //   {
    //     if (physVol.children().size() == 0)
    //       writeLVToPV(pvName, lvName, pvID, lvID);
    //   }
    // }
    // else
    {
        // writeSolidToLV(lvName, solidName, (*m_Visible)[lvIndex], matID, lvID+lvCount, solidID);
        // writeLVToPV(pvName, lvName, pvID+pvCount, lvID+lvCount);
        // writeSolidToPV(pvName, solidName, (*m_Visible)[lvIndex], matID, pvID + pvCount, solidID);
    }(*m_LVCount)[lvIndex]++;
    (*m_PVCount)[pvIndex]++;
  }
}

void dd4hep2FBXWriter::addModels(TGeoNode *node, int replica, unsigned long long pvIndex)
{
  unsigned long long pvID = (*m_PVID)[pvIndex];
  // unsigned int pvCount = (*m_PVCount)[pvIndex];  // 空的 0
  std::string pvName = m_placedVolName[pvIndex];
  int lvIndex = pvIndex; // std::find(lvStore->begin(), lvStore->end(), logVol) - lvStore->begin();
  unsigned long long lvID = (*m_LVID)[lvIndex];
  // unsigned int lvCount = (*m_LVCount)[lvIndex];  // 空的 0
  std::string lvName = m_volName[lvIndex];
  // if ((*m_LVUnique)[lvIndex]) writeMaterialNode(lvIndex, (*m_PVName)[pvIndex]);

  // if (m_UsePrototypes)
  // {
  //   // if (lvCount == 0)
  //   {
  //     // if (!(*m_LVUnique)[lvIndex])
  writeLVModelNode(node->GetVolume(), lvName, lvID);
  //   }
  // if (pvCount == 0)
  writePVModelNode(node, pvName, pvID);
  // }
  // else
  // {
  //   // DIVOT writeLVModelNode((*lvStore)[lvIndex], lvName, lvID+lvCount);
  //   // writePVModelNode(physVol, pvName, pvID + pvCount);
  // }
  // (*m_LVCount)[lvIndex]++;
  // (*m_PVCount)[pvIndex]++;
  // }
}

void dd4hep2FBXWriter::writeLVModelNode(TGeoVolume *vol, const std::string lvName, unsigned long long lvID)
{
  m_File << "\t; LogVol " << vol->GetName() << " with solid " << vol->GetShape()->GetName() << std::endl
         << "\tModel: " << lvID << ", \"Model::lv_" << lvName << "\", \"Null\" {" << std::endl
         << "\t\tVersion: 232" << std::endl
         << "\t\tProperties70:  {" << std::endl
         << "\t\t}" << std::endl
         << "\t\tShading: T" << std::endl
         << "\t\tCulling: \"CullingOff\"" << std::endl
         << "\t}" << std::endl;
}

void dd4hep2FBXWriter::writePVModelNode(TGeoNode *node, const std::string pvName, unsigned long long pvID)
{
  TGeoMatrix *matrix = node->GetMatrix();
  const Double_t *tv = matrix->GetTranslation();

  // G4RotationMatrix* rot = physVol->GetObjectRotation();
  // G4ThreeVector move = physVol->GetObjectTranslation();
  // FBX uses the Tait-Bryan version of the Euler angles (X then Y then Z rotation)
  TGeoRotation *mat = (TGeoRotation *)node->GetMatrix();
  double phi = 0, theta = 0, psi = 0;
  mat->GetAngles(phi, theta, psi); // all angles in degrees

  /*// Construct from three Euler angles (in radians).
  CLHEP::HepRotation *rot = new
      CLHEP::HepRotation(phi*M_PI/180.0, theta*M_PI/180.0, psi*M_PI/180.0);

  double yaw = std::atan2(rot->yx(), rot->xx()) * 180.0 / M_PI;
  if (fabs(yaw) < 1.0E-12) yaw = 0.0;
  double pitch = -std::asin(rot->zx()) * 180.0 / M_PI;
  if (fabs(pitch) < 1.0E-12) pitch = 0.0;
  double roll = std::atan2(rot->zy(), rot->zz()) * 180.0 / M_PI;
  if (fabs(roll) < 1.0E-12) roll = 0.0;
  */
  m_File << "\t; PhysVol " << node->GetName();
  // if (physVol->IsReplicated()) {
  //   m_File << " (replicated: copy " << physVol->GetCopyNo() << ")";
  // }
  m_File << ", placing LogVol " << node->GetVolume()->GetName() << std::endl
         << "\tModel: " << pvID << ", \"Model::" << pvName << "\", \"Null\" {" << std::endl
         << "\t\tVersion: 232" << std::endl
         << "\t\tProperties70:  {" << std::endl
         << "\t\t\tP: \"Lcl Translation\", \"Lcl Translation\", \"\", \"A\"," << tv[0] << "," << tv[1] << "," << tv[2] << std::endl
         << "\t\t\tP: \"Lcl Rotation\", \"Lcl Rotation\", \"\", \"A\"," << phi << "," << theta << "," << psi << std::endl
         << "\t\t}" << std::endl
         << "\t\tShading: T" << std::endl
         << "\t\tCulling: \"CullingOff\"" << std::endl
         << "\t}" << std::endl;
}

void dd4hep2FBXWriter::writeMaterialNode(int lvIndex, const std::string matName)
{
  // G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  // Volume logVol = m_vol[lvIndex];
  unsigned long long matID = (*m_MatID)[lvIndex];
  float alpha = m_color[i]->GetAlpha(),
        red = m_color[i]->GetRed(),
        green = m_color[i]->GetGreen(),
        blue = m_color[i]->GetBlue();

  bool visible = m_visible[lvIndex];
  string materialName = logVol.material().name();
  // Hide volumes that contain vacuum, air or gas
  // if (materialName == "Vacuum")
  //   visible = false;
  // if (materialName == "G4_AIR")
  //   visible = false;
  // Hide volumes that are invisible in the GEANT4 geometry
  m_File << "\t; Color for LogVol " << m_volName[lvIndex] << std::endl
         << "\tMaterial: " << matID << ", \"Material::" << matName << "\", \"\" {" << std::endl
         << "\t\tVersion: 102" << std::endl
         << "\t\tProperties70:  {" << std::endl
         << "\t\t\tP: \"ShadingModel\", \"KString\", \"\", \"\", \"phong\"" << std::endl
         << "\t\t\tP: \"DiffuseColor\", \"RGBColor\", \"Color\", \"A\"," << red << "," << green << "," << blue << std::endl
         << "\t\t\tP: \"TransparentColor\", \"RGBColor\", \"Color\", \"A\",1,1,1" << std::endl
         << "\t\t\tP: \"TransparencyFactor\", \"double\", \"Number\", \"\"," << (visible ? 1.0 - alpha : 1) << std::endl
         << "\t\t}" << std::endl
         << "\t}" << std::endl;
}

double GetCutZ(const CLHEP::Hep3Vector &p, const CLHEP::Hep3Vector fLowNorm, const CLHEP::Hep3Vector fHighNorm, double fDz)
{
  double newz = p.z(); // p.z() should be either +fDz or -fDz
  if (p.z() < 0)
  {
    if (fLowNorm.z() != 0.)
    {
      newz = -fDz - (p.x() * fLowNorm.x() + p.y() * fLowNorm.y()) / fLowNorm.z();
    }
  }
  else
  {
    if (fHighNorm.z() != 0.)
    {
      newz = fDz - (p.x() * fHighNorm.x() + p.y() * fHighNorm.y()) / fHighNorm.z();
    }
  }
  return newz;
}
GPolyhedron *getPolyhedron(Solid solid)
{

  // https://root.cern.ch/doc/master/classTGDMLWrite.html
  string solidtype = solid.type();
  if (solidtype == "TGeoShapeAssembly") // ShapelessSolid
  {
  }
  else if (solidtype == "TGeoScaledShape") // solid Scale
  {
    Scale tmp(solid);
    double scale_x = tmp.scale_x();
    // return new GPolyhedronCons(tmp.rMin1(), tmp.rMax1(), tmp.rMin2(), tmp.rMax2(), tmp.dZ(), 0, twopi);
  }
  else if (solidtype == "TGeoBBox")
  {
    return new GPolyhedronBox(Box(solid).x(), Box(solid).y(), Box(solid).z());
  }
  else if (solidtype == "TGeoHalfSpace") // HalfSpace
  {
  }
  else if (solidtype == "TGeoCone") // cone
  {
    Cone tmp(solid);
    return new GPolyhedronCons(tmp.rMin1(), tmp.rMax1(), tmp.rMin2(), tmp.rMax2(), tmp.dZ(), 0, twopi);
  }
  else if (solidtype == "TGeoPcon") //  Polycone
  {
    // https://apc.u-paris.fr/~franco/g4doxy/html/G4Polycone_8cc-source.html#l00898
    Polycone tmp(solid);
    const int nz = tmp.zPlaneZ().size();
    double Z_values[nz];
    double Rmin[nz];
    double Rmax[nz];
    for (int i = 0; i < nz; i++)
    {
      Z_values[i] = tmp.z(i);
      Rmin[i] = tmp.rMin(i);
      Rmax[i] = tmp.rMax(i);
    }
    return new GPolyhedronPcon(tmp.startPhi(), // original_parameters->Start_angle,
                               tmp.deltaPhi(), // original_parameters->Opening_angle,
                               nz,             // original_parameters->Num_z_planes,
                               Z_values,       // original_parameters->Z_values,
                               Rmin,           // original_parameters->Rmin,
                               Rmax);          // original_parameters->Rmax );
  }
  else if (solidtype == "TGeoConeSeg") // ConeSegment
  {
    ConeSegment tmp(solid);
    return new GPolyhedronCons(tmp.rMin1(), tmp.rMax1(), tmp.rMin2(), tmp.rMax2(), tmp.dZ(), tmp.startPhi(), tmp.endPhi());
  }
  else if (solidtype == "TGeoTubeSeg") // Tube___________
  {
    Tube tmp(solid);
    return new GPolyhedronTubs(tmp.rMin(), tmp.rMax(), tmp.dZ(), tmp.startPhi(), tmp.endPhi());
  }
  else if (solidtype == "TGeoCtub") // CutTube
  {
    // https://apc.u-paris.fr/~franco/g4doxy/html/G4CutTubs_8cc-source.html#l01975
    CutTube tmp(solid);
    typedef double double3[3];
    typedef int int4[4];
    double fRMin = tmp.rMin(), fRMax = tmp.rMax(), fDz = tmp.dZ(), fSPhi = tmp.startPhi(), fDPhi = tmp.endPhi();
    double kCarTolerance = 0; // Cached geometrical tolerance ??

    GPolyhedron *ph = new GPolyhedron;
    GPolyhedron *ph1 = new GPolyhedronTubs(fRMin, fRMax, fDz, fSPhi, fDPhi); // G4Tubs::CreatePolyhedron();
    int nn = ph1->GetNoVertices();
    int nf = ph1->GetNoFacets();
    double3 *xyz = new double3[nn]; // number of nodes
    int4 *faces = new int4[nf];     // number of faces

    std::vector<double> lowNormal = tmp.lowNormal();
    std::vector<double> highNormal = tmp.highNormal();
    const CLHEP::Hep3Vector fLowNorm(lowNormal[0], lowNormal[1], lowNormal[2]);
    const CLHEP::Hep3Vector fHighNorm(highNormal[0], highNormal[1], highNormal[2]);
    for (int i = 0; i < nn; i++)
    {
      xyz[i][0] = ph1->GetVertex(i + 1).x();
      xyz[i][1] = ph1->GetVertex(i + 1).y();
      double tmpZ = ph1->GetVertex(i + 1).z();
      if (tmpZ >= fDz - kCarTolerance)
      {
        xyz[i][2] = GetCutZ(CLHEP::Hep3Vector(xyz[i][0], xyz[i][1], fDz), fLowNorm, fHighNorm, fDz);
      }
      else if (tmpZ <= -fDz + kCarTolerance)
      {
        xyz[i][2] = GetCutZ(CLHEP::Hep3Vector(xyz[i][0], xyz[i][1], -fDz), fLowNorm, fHighNorm, fDz);
      }
      else
      {
        xyz[i][2] = tmpZ;
      }
    }
    int iNodes[4];
    int *iEdge = 0;
    int n;
    for (int i = 0; i < nf; i++)
    {
      ph1->GetFacet(i + 1, n, iNodes, iEdge);
      for (int k = 0; k < n; k++)
      {
        faces[i][k] = iNodes[k];
      }
      for (int k = n; k < 4; k++)
      {
        faces[i][k] = 0;
      }
    }
    ph->createPolyhedron(nn, nf, xyz, faces);

    delete[] xyz;
    delete[] faces;
    delete ph1;

    return ph;
  }
  else if (solidtype == "TGeoEltu") // EllipticalTube
  {
    EllipticalTube tmp(solid);
    GPolyhedronTube *eTube = new GPolyhedronTube(0., 1., tmp.dZ());
    // apply non-uniform scaling...
    // undefined reference to `HepGeom::operator*(HepGeom::Transform3D const&, HepGeom::Vector3D<double> const&)'
    // eTube->Transform(HepGeom::Scale3D(tmp.a(), tmp.b(),1.));
    // return  eTube;
  }
  else if (solidtype == "TGeoTubeSeg") // TwistedTube__________
  {
    // Trap tmp(solid);
    // return new GPolyhedronTubs(tmp.rMin(), tmp.rMax(), tmp.dZ(), tmp.startPhi(), tmp.endPhi()) ;
  }
  else if (solidtype == "TGeoTrap") // Trap
  {
    Trap tmp(solid);
    return new GPolyhedronTrap(tmp.dZ(), tmp.theta(), tmp.phi(),
                               tmp.high1(), tmp.bottomLow1(), tmp.topLow1(), tmp.alpha1(),
                               tmp.high2(), tmp.bottomLow2(), tmp.topLow2(), tmp.alpha2());
  }
  else if (solidtype == "TGeoTrd1") // Trd1
  {
    Trd1 tmp(solid);
    return new GPolyhedronTrd1(tmp.dX1(), tmp.dX2(), tmp.dY(), tmp.dZ());
  }
  else if (solidtype == "TGeoTrd2") // Trd2
  {
    Trd2 tmp(solid);
    return new GPolyhedronTrd2(tmp.dX1(), tmp.dX2(), tmp.dY1(), tmp.dY2(), tmp.dZ());
  }
  else if (solidtype == "TGeoTorus") // Torus
  {
    Torus tmp(solid);
    return new GPolyhedronTorus(tmp.rMin(), tmp.rMax(), tmp.r(), tmp.startPhi(), tmp.deltaPhi());
  }
  else if (solidtype == "TGeoSphere") // Sphere
  {
    Sphere tmp(solid);
    return new GPolyhedronSphere(tmp.rMin(), tmp.rMax(), tmp.startPhi(), tmp.endPhi(), tmp.startTheta(), tmp.endTheta());
  }
  else if (solidtype == "TGeoParaboloid") // Paraboloid
  {
    Paraboloid tmp(solid);
    return new GPolyhedronParaboloid(tmp.rLow(), tmp.rHigh(), tmp.dZ(), 0., twopi);
  }
  else if (solidtype == "TGeoHype") // Hyperboloid
  {
    Hyperboloid tmp(solid);
    return new GPolyhedronHype(tmp.rMin(), tmp.rMax(),
                               tmp.stereoInner(), tmp.stereoOuter(), tmp.dZ());
  }
  else if (solidtype == "TGeoPgon") // PolyhedraRegular______
  {
  }
  else if (solidtype == "TGeoPgon") // Polyhedra_____
  {
  }
  else if (solidtype == "TGeoXtru") // ExtrudedPolygon
  {
  }
  else if (solidtype == "TGeoArb8") // EightPointSolid
  {
  }
  else if (solidtype == "TGeoTessellated") // TessellatedSolid
  {
    // https://apc.u-paris.fr/~franco/g4doxy/html/classG4TessellatedSolid.html#8f487866dbc90c2fddf520fcc7289359
  }

  return nullptr;
}

void dd4hep2FBXWriter::writeGeometryNode(Solid solid, const std::string solidName, unsigned long long solidID)
{
  std::cout << __LINE__ << " solid.type:" << solid.type() << /*" IsComposite:" << solid.IsComposite() << */ std::endl;
  // IntersectionSolid, SubtractionSolid, UnionSolid are retrived from BooleanSolid
  // all their types are TGeoCompositeShape
  string solidtype = solid.type();
  if (solidtype == "TGeoCompositeShape")
  {
    HepPolyhedron *polyhedron = getBooleanSolidPolyhedron(solid);
    GPolyhedron *polyh = new GPolyhedron(*polyhedron);
    writePolyhedron(solid, polyh, solidName, solidID);
    delete polyhedron;
    delete polyh;
  }
  else
  {
    // auto a=Polyhedra(solid);
    std::cout << __LINE__ << " name:" << solid.name() << " type:" << solid.type() << std::endl;
    std::cout << __LINE__ << "solid info: " << solid.toString() << std::endl;
    writePolyhedron(solid, getPolyhedron(solid), solidName, solidID);
  }
}

HepPolyhedron *dd4hep2FBXWriter::getBooleanSolidPolyhedron(Solid solid)
{

  BooleanSolid boSolid = (BooleanSolid)solid;
  Solid solidA = boSolid.rightShape();
  Solid solidB = boSolid.leftShape();

  HepPolyhedron *polyhedronA = NULL;
  string solidAtype = solidA.type();
  if ((solidAtype == "TGeoCompositeShape"))
  {
    polyhedronA = getBooleanSolidPolyhedron(solidA);
  }
  else
  {
    polyhedronA = new HepPolyhedron(*(getPolyhedron(solidA)));
  }

  HepPolyhedron *polyhedronB = NULL;
  string solidBtype = solidB.type();
  if ((solidBtype == "TGeoCompositeShape"))
  {
    polyhedronB = getBooleanSolidPolyhedron(solidB);
  }
  else
  {
    polyhedronB = new HepPolyhedron(*(getPolyhedron(solidB)));
  }

  /// home/wln/DD4hep_source/DDCore/src/ShapeUtilities.cpp
  TGeoCompositeShape *sh = (TGeoCompositeShape *)&(*solid);
  TGeoBoolNode *boolean = sh->GetBoolNode();
  TGeoBoolNode::EGeoBoolType oper = boolean->GetBooleanOperator();

  HepPolyhedron *result = new HepPolyhedron();
  if (oper == TGeoBoolNode::kGeoSubtraction)
    *result = polyhedronA->subtract(*polyhedronB);
  else if (oper == TGeoBoolNode::kGeoUnion)
    *result = polyhedronA->add(*polyhedronB);
  else if (oper == TGeoBoolNode::kGeoIntersection)
    *result = polyhedronA->intersect(*polyhedronB);
  else
  {
    std::cerr << "getBooleanSolidPolyhedron(): Unrecognized boolean solid " << solid.name() << " of type " << solid.type() << std::endl;
  }
  delete polyhedronA;
  delete polyhedronB;
  return result;
}

void dd4hep2FBXWriter::writePolyhedron(Solid solid, GPolyhedron *polyhedron, const std::string name,
                                       unsigned long long solidID)
{
  if (polyhedron)
  {
    polyhedron->SetNumberOfRotationSteps(120);
    m_File << "\t; Solid " << solid.name() << " of type " << solid.type() << std::endl
           << "\tGeometry: " << solidID << ", \"Geometry::" << name << "\", \"Mesh\" {" << std::endl
           << "\t\tVertices: *" << polyhedron->GetNoVertices() * 3 << " {" << std::endl
           << "\t\t\ta: ";
    std::streampos startOfLine = m_File.tellp();
    for (int j = 1; j <= polyhedron->GetNoVertices(); ++j)
    {
      m_File << (j == 1 ? "" : ",") << polyhedron->GetVertex(j).x() << "," << polyhedron->GetVertex(j).y() << "," << polyhedron->GetVertex(j).z();
      if (m_File.tellp() - startOfLine > 100)
      {
        startOfLine = m_File.tellp();
        m_File << std::endl
               << "\t\t\t\t";
      }
    }
    m_File << std::endl
           << "\t\t}" << std::endl;

    std::vector<int> vertices;
    for (int k = 1; k <= polyhedron->GetNoFacets(); ++k)
    {
      bool notLastEdge = true;
      int ndx = -1, edgeFlag = 1;
      do
      {
        notLastEdge = polyhedron->GetNextVertexIndex(ndx, edgeFlag);
        if (notLastEdge)
        {
          vertices.push_back(ndx - 1);
        }
        else
        {
          vertices.push_back(-ndx);
        }
      } while (notLastEdge);
    }
    m_File << "\t\tPolygonVertexIndex: *" << vertices.size() << " {" << std::endl
           << "\t\t\ta: ";
    startOfLine = m_File.tellp();
    for (unsigned int j = 0; j < vertices.size(); ++j)
    {
      m_File << (j == 0 ? "" : ",") << vertices[j];
      if (m_File.tellp() - startOfLine > 100)
      {
        startOfLine = m_File.tellp();
        m_File << std::endl
               << "\t\t\t\t";
      }
    }
    m_File << std::endl
           << "\t\t}" << std::endl;

    m_File << "\t\tGeometryVersion: 124" << std::endl
           << "\t\tLayerElementNormal: 0 {" << std::endl
           << "\t\t\tVersion: 101" << std::endl
           <<
        // "\t\t\tName: \"\"" << std::endl <<
        "\t\t\tMappingInformationType: \"ByPolygonVertex\"" << std::endl
           << "\t\t\tReferenceInformationType: \"Direct\"" << std::endl
           << "\t\t\tNormals: *" << vertices.size() * 3 << " {" << std::endl
           << "\t\t\t\ta: ";
    startOfLine = m_File.tellp();
    unsigned int j = 0;
    for (int k = 1; k <= polyhedron->GetNoFacets(); ++k)
    {
      HepGeom::Normal3D<double> normal = polyhedron->GetUnitNormal(k);
      do
      {
        m_File << (j == 0 ? "" : ",") << normal.x() << "," << normal.y() << "," << normal.z();
        if (m_File.tellp() - startOfLine > 100)
        {
          startOfLine = m_File.tellp();
          m_File << std::endl
                 << "\t\t\t\t";
        }
      } while (vertices[j++] >= 0);
    }
    m_File << std::endl
           << "\t\t\t}" << std::endl
           << "\t\t}" << std::endl
           << "\t\tLayerElementMaterial: 0 {" << std::endl
           << "\t\t\tVersion: 101" << std::endl
           <<
        // "\t\t\tName: \"\"" << std::endl <<
        "\t\t\tMappingInformationType: \"AllSame\"" << std::endl
           << "\t\t\tReferenceInformationType: \"IndexToDirect\"" << std::endl
           << "\t\t\tMaterials: *1 {" << std::endl
           << "\t\t\t\ta: 0" << std::endl
           << "\t\t\t}" << std::endl
           << "\t\t}" << std::endl
           << "\t\tLayer: 0 {" << std::endl
           << "\t\t\tVersion: 100" << std::endl
           << "\t\t\tLayerElement:  {" << std::endl
           << "\t\t\t\tType: \"LayerElementNormal\"" << std::endl
           << "\t\t\t\tTypedIndex: 0" << std::endl
           << "\t\t\t}" << std::endl
           << "\t\t\tLayerElement:  {" << std::endl
           << "\t\t\t\tType: \"LayerElementMaterial\"" << std::endl
           << "\t\t\t\tTypedIndex: 0" << std::endl
           << "\t\t\t}" << std::endl
           << "\t\t}" << std::endl
           << "\t}" << std::endl;
  }
  else
  {
    std::cerr << "Polyhedron representation of solid " << name << " cannot be created" << std::endl;
  }
}

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

void dd4hep2FBXWriter::writeSolidToLV(const std::string lvName, const std::string solidName, bool visible,
                                      unsigned long long matID, unsigned long long lvID, unsigned long long solidID)
{
  m_File << "\t; Solid Geometry::" << solidName << ", LogVol Model::lv_" << lvName << std::endl
         << "\t" << (visible ? "" : "; ") << "C: \"OO\"," << solidID << "," << lvID << std::endl
         << std::endl
         << "\t; Color Material::" << lvName << ", LogVol Model::lv_" << lvName << std::endl
         << "\t" << (visible ? "" : "; ") << "C: \"OO\"," << matID << "," << lvID << std::endl
         << std::endl;
}

void dd4hep2FBXWriter::writeSolidToPV(const std::string pvName, const std::string solidName, bool visible,
                                      unsigned long long matID, unsigned long long pvID, unsigned long long solidID)
{
  m_File << "\t; Solid Geometry::" << solidName << ", PhysVol Model::" << pvName << std::endl
         << "\t" << (visible ? "" : "; ") << "C: \"OO\"," << solidID << "," << pvID << std::endl
         << std::endl
         << "\t; Color Material::" << pvName << ", PhysVol Model::" << pvName << std::endl
         << "\t" << (visible ? "" : "; ") << "C: \"OO\"," << matID << "," << pvID << std::endl
         << std::endl;
}

void dd4hep2FBXWriter::writeLVToPV(const std::string pvName, const std::string lvName, unsigned long long pvID,
                                   unsigned long long lvID)
{
  m_File << "\t; LogVol Model::lv_" << lvName << ", PhysVol Model::" << pvName << std::endl
         << "\tC: \"OO\"," << lvID << "," << pvID << std::endl
         << std::endl;
}

void dd4hep2FBXWriter::writePVToParentLV(const std::string pvNameDaughter, const std::string lvName,
                                         unsigned long long pvIDDaughter, unsigned long long lvID)
{
  m_File << "\t; PhysVol Model::" << pvNameDaughter << ", parent LogVol Model::lv_" << lvName << std::endl
         << "\tC: \"OO\"," << pvIDDaughter << "," << lvID << std::endl
         << std::endl;
}

void dd4hep2FBXWriter::writePVToParentPV(const std::string pvNameDaughter, const std::string pvName,
                                         unsigned long long pvIDDaughter, unsigned long long pvID)
{
  m_File << "\t; PhysVol Model::" << pvNameDaughter << ", parent PhysVol Model::" << pvName << std::endl
         << "\tC: \"OO\"," << pvIDDaughter << "," << pvID << std::endl
         << std::endl;
}