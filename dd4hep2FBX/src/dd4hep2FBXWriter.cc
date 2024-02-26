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

void dd4hep2FBXWriter::getAllChildren(DetElement det)
{
  m_childrenDet.push_back(det);
  m_childrenVol.push_back(det.volume());
  m_childrenSolid.push_back(det.solid());

  m_DetName.push_back(det.name());
  m_VolName.push_back(det.volume().name());
  // det.volume()->setName("adf");
  m_SolidName.push_back(det.solid().name());
  // det.solid().setName("adf");
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
  m_PVID = new std::vector<unsigned long long>(m_childrenDet.size(), 0x0000010000000000LL);
  m_LVID = new std::vector<unsigned long long>(m_childrenVol.size(), 0x000000C000000000LL);
  m_SolidID = new std::vector<unsigned long long>(m_childrenSolid.size(), 0x0000008000000000LL);
  m_MatID = new std::vector<unsigned long long>(m_childrenVol.size(), 0x0000004000000000LL);
  m_Visible = new std::vector<bool>(m_childrenVol.size(), false);

  m_File << "Objects:  {" << std::endl;
  for (unsigned int solidIndex = 0; solidIndex < m_childrenSolid.size(); ++solidIndex)
  {
    (*m_SolidID)[solidIndex] += 0x0000000001000000LL * solidIndex;
    if (m_SolidName[solidIndex].length() > 0)
    {
      std::cout << __LINE__ << " " << m_childrenSolid.size() << " " << m_SolidID->size()
                << m_SolidName[solidIndex] << " "
                << " index " << solidIndex << " "
                << (*m_SolidID)[solidIndex] << std::endl;
      writeGeometryNode(m_childrenSolid[solidIndex], m_SolidName[solidIndex], (*m_SolidID)[solidIndex]);
    }
  }

  // write materials
  for (unsigned int lvIndex = 0; lvIndex < m_childrenVol.size(); ++lvIndex)
  {
    (*m_MatID)[lvIndex] += 0x0000000001000000LL * lvIndex;
    (*m_LVID)[lvIndex] += 0x0000000001000000LL * lvIndex;
    // if (m_VolName[lvIndex].length() > 0) {
    // if (!(*m_LVUnique)[lvIndex])
    std::cout << __LINE__ << " materials_childVol_index: " << lvIndex << std::endl;
    writeMaterialNode(lvIndex, m_VolName[lvIndex]);

    // }
  }

  for (unsigned int pvIndex = 0; pvIndex < m_childrenDet.size(); ++pvIndex)
  {
    (*m_PVID)[pvIndex] += 0x0000000001000000LL * pvIndex;
  }

  //

  // m_PVCount->assign(m_childrenDet.size(), 0);
  // std::cout << __LINE__ << " addmodels " << std::endl;
  // m_LVCount->assign(m_childrenVol.size(), 0);
  // m_SolidCount->assign(m_childrenSolid.size(), 0);
  for (unsigned int i = 0; i < m_childrenDet.size(); i++)
  {
    // if(m_childrenDet[i] == m_world)
    // std::cout << __LINE__ << " addmodels " << i << std::endl;
    addModels(m_childrenDet[i], 0, i);
    // else
    //   addModels(m_childrenDet[i], 1);
  }
  m_File << "}" << std::endl
         << std::endl;

  // Recursively write the connections among the solid and logical/physical volume elements
  m_File << "Connections:  {" << std::endl;
  addConnections(m_world, 0);

  // int pvIndex = std::find(m_childrenDet.begin(), m_childrenDet.end(), m_world) - m_childrenDet.begin();
  m_File << "\t; Physical volume Model::" << m_world.name() << " to Model::RootNode" << std::endl
         << "\tC: \"OO\"," << 0 << ",0" << std::endl
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
  int pvIndex = std::find(m_childrenDet.begin(), m_childrenDet.end(), physVol) - m_childrenDet.begin();
  ; //
  /// std::find(pvStore->begin(), pvStore->end(), physVol) - pvStore->begin();
  unsigned long long pvID = (*m_PVID)[pvIndex];
  // unsigned int pvCount = (*m_PVCount)[pvIndex];
  std::string pvName = m_DetName[pvIndex];
  int lvIndex = pvIndex; // std::find(lvStore->begin(), lvStore->end(), logVol) - lvStore->begin();
  unsigned long long lvID = (*m_LVID)[lvIndex];
  // unsigned int lvCount = (*m_LVCount)[lvIndex];
  std::string lvName = m_VolName[lvIndex];
  Children mchildren = physVol.children();
  // for (long unsigned int daughter = 0; daughter < mchildren.size(); ++daughter) {
  for (const auto &[name, physVolDaughter] : mchildren)
  {
    // DetElement physVolDaughter = mchildren.at(daughter).second ;// ->GetDaughter(daughter);
    int pvIndexDaughter = std::find(m_childrenDet.begin(), m_childrenDet.end(), physVolDaughter) - m_childrenDet.begin();
    unsigned long long pvIDDaughter = (*m_PVID)[pvIndexDaughter];
    // unsigned int pvCountDaughter = (*m_PVCount)[pvIndexDaughter];
    // for (int j = 0; j < physVolDaughter->GetMultiplicity(); ++j) {
    //   if (m_UsePrototypes) {
    //     if ((replica == 0) && (j == 0) && (lvCount == 0) && (pvCountDaughter == 0)) {
    //       if ((*m_LVUnique)[lvIndex]) {
    writePVToParentPV(m_DetName[pvIndexDaughter], pvName, pvIDDaughter, pvID);
    //     } else {
    //       writePVToParentLV(m_DetName[pvIndexDaughter], lvName, pvIDDaughter, lvID);
    //     }
    //   }
    // } else {
    //   //writePVToParentLV(m_DetName[pvIndexDaughter], lvName, pvIDDaughter+0x00010000*j+pvCountDaughter, lvID+0x00010000*replica+lvCount);
    //   writePVToParentPV(m_DetName[pvIndexDaughter], pvName, pvIDDaughter + 0x00010000 * j + pvCountDaughter,
    //                     pvID + 0x00010000 * replica + pvCount);
    // }
    addConnections(physVolDaughter, 0);
    // }
  }

  // Write the Geometry-LogVolModel, Material-LogVolModel and PhysVolModel-LogVolModel
  // connections as we ascend the recursive tree
  Solid solid = logVol.solid();
  int solidIndex = pvIndex; // std::find(solidStore->begin(), solidStore->end(), solid) - solidStore->begin();
  unsigned long long solidID = (*m_SolidID)[solidIndex];
  unsigned long long matID = (*m_MatID)[lvIndex];
  std::string solidName = m_SolidName[solidIndex];
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
    {
      // if (lvCount == 0)
      // {
      // if ((*m_LVUnique)[lvIndex])
      { // bypass the singleton logical volume
        // writeSolidToPV(pvName, solidName, (*m_Visible)[lvIndex], matID, pvID, solidID);
      }
      // else
      // {
      writeSolidToLV(lvName, solidName, (*m_Visible)[lvIndex], matID, lvID, solidID);
      // }
      // }
      // if (pvCount == 0)
      {
        // if (!(*m_LVUnique)[lvIndex])
        writeLVToPV(pvName, lvName, pvID, lvID);
      }
    }
    // else
    {
      // writeSolidToLV(lvName, solidName, (*m_Visible)[lvIndex], matID, lvID+lvCount, solidID);
      // writeLVToPV(pvName, lvName, pvID+pvCount, lvID+lvCount);
      // writeSolidToPV(pvName, solidName, (*m_Visible)[lvIndex], matID, pvID + pvCount, solidID);
    }
    // (*m_LVCount)[lvIndex]++;
    // (*m_PVCount)[pvIndex]++;
  }
}

void dd4hep2FBXWriter::addModels(DetElement physVol, int replica, unsigned long long pvIndex)
{
  // Descend to the leaves of the tree
  // DetElement logVol = physVol.volume();
  // for (int daughter = 0; daughter < logVol->GetNoDaughters(); ++daughter) {
  //   G4VPhysicalVolume* physVolDaughter = logVol->GetDaughter(daughter);
  //   for (int j = 0; j < physVolDaughter->GetMultiplicity(); ++j) {
  //     addModels(physVolDaughter, j);
  //   }
  // }

  // Write the physical- and logical-volume models as we ascend the recursive tree
  // G4PhysicalVolumeStore* pvStore = G4PhysicalVolumeStore::GetInstance();
  // G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  // int pvIndex = std::find(pvStore->begin(), pvStore->end(), physVol) - pvStore->begin();
  unsigned long long pvID = (*m_PVID)[pvIndex];
  // unsigned int pvCount = (*m_PVCount)[pvIndex];  // 空的 0
  std::string pvName = m_DetName[pvIndex];
  int lvIndex = pvIndex; // std::find(lvStore->begin(), lvStore->end(), logVol) - lvStore->begin();
  unsigned long long lvID = (*m_LVID)[lvIndex];
  // unsigned int lvCount = (*m_LVCount)[lvIndex];  // 空的 0
  std::string lvName = m_VolName[lvIndex];
  // if ((*m_LVUnique)[lvIndex]) writeMaterialNode(lvIndex, (*m_PVName)[pvIndex]);
  /*if (physVol->IsReplicated()) {
    G4VSolid* solid = logVol->GetSolid();
    G4SolidStore* solidStore = G4SolidStore::GetInstance();
    int solidIndex = std::find(solidStore->begin(), solidStore->end(), solid) - solidStore->begin();
    unsigned long long solidID = (*m_SolidID)[solidIndex];
    EAxis axis;
    G4int nReplicas;
    G4double width;
    G4double offset;
    G4bool consuming;
    physVol->GetReplicationData(axis, nReplicas, width, offset, consuming);
    physVol->SetCopyNo(replica);
    G4VPVParameterisation* physParameterisation = physVol->GetParameterisation();
    if (physParameterisation) { // parameterised volume
      G4VSolid* solidReplica = physParameterisation->ComputeSolid(replica, physVol);
      physParameterisation->ComputeTransformation(replica, physVol);
      solidReplica->ComputeDimensions(physParameterisation, replica, physVol);
      if (!(*solidReplica == *solid)) {
        std::string solidName = (*m_SolidName)[solidIndex];
        solidName.append("_R");
        solidName.append(std::to_string(replica));
        writeGeometryNode(solidReplica, solidName, solidID + 0x00010000 * replica);
      }
      if (m_UsePrototypes && (*solidReplica == *solid)) {
        if ((replica == 0) && (lvCount == 0)) {
          if (!(*m_LVUnique)[lvIndex]) writeLVModelNode((*lvStore)[lvIndex], lvName, lvID);
        }
      } else {
        // DIVOT lvName.append("_R");
        // DIVOT lvName.append(std::to_string(replica));
        // DIVOT writeLVModelNode((*lvStore)[lvIndex], lvName, lvID+0x00010000*replica+lvCount);
      }
      pvName.append("_R");
      pvName.append(std::to_string(replica));
      writePVModelNode(physVol, pvName, pvID + 0x00010000 * replica + pvCount);
    } else { // plain replicated volume
      G4RotationMatrix* originalRotation = physVol->GetRotation();
      G4ThreeVector translation; // No translation
      G4RotationMatrix rotation; // No rotation
      switch (axis) {
        default:
        case kXAxis:
          translation.setX(width * (replica - 0.5 * (nReplicas - 1)));
          physVol->SetTranslation(translation);
          break;
        case kYAxis:
          translation.setY(width * (replica - 0.5 * (nReplicas - 1)));
          physVol->SetTranslation(translation);
          break;
        case kZAxis:
          translation.setZ(width * (replica - 0.5 * (nReplicas - 1)));
          physVol->SetTranslation(translation);
          break;
        case kRho:
          if (solid->GetEntityType() == "G4Tubs") {
            double originalRMin = ((G4Tubs*)solid)->GetInnerRadius();
            double originalRMax = ((G4Tubs*)solid)->GetOuterRadius();
            ((G4Tubs*)solid)->SetInnerRadius(offset + width * replica);
            ((G4Tubs*)solid)->SetOuterRadius(offset + width * (replica + 1));
            std::string solidName = (*m_SolidName)[solidIndex];
            solidName.append("_R");
            solidName.append(std::to_string(replica));
            writeGeometryNode(solid, solidName, (*m_SolidID)[solidIndex] + 0x00010000 * replica);
            ((G4Tubs*)solid)->SetInnerRadius(originalRMin);
            ((G4Tubs*)solid)->SetOuterRadius(originalRMax);
          } else if (replica == 0) {
            std::err << "Built-in volumes replicated along radius for " << solid->GetEntityType() <<
                      " (solid " << solid->GetName() << ") are not visualisable." << std::endl;
          }
          break;
        case kPhi:
          physVol->SetRotation(&(rotation.rotateZ(-(offset + (replica + 0.5) * width))));
          break;
      }
      if (m_UsePrototypes && !((axis == kRho) && (solid->GetEntityType() == "G4Tubs"))) {
        if ((replica == 0) && (lvCount == 0)) {
          if (!(*m_LVUnique)[lvIndex]) writeLVModelNode((*lvStore)[lvIndex], lvName, lvID);
        }
      } else {
        // DIVOT lvName.append("_R");
        // DIVOT lvName.append(std::to_string(replica));
        // DIVOT writeLVModelNode((*lvStore)[lvIndex], lvName, lvID+0x00010000*replica+lvCount);
      }
      pvName.append("_R");
      pvName.append(std::to_string(replica));
      writePVModelNode(physVol, pvName, pvID + 0x00010000 * replica + pvCount);
      if (axis == kPhi) physVol->SetRotation(originalRotation);
    }
  }
  else*/
  // {
  // if (m_UsePrototypes)
  // {
  //   // if (lvCount == 0)
  //   {
  //     // if (!(*m_LVUnique)[lvIndex])
  //       writeLVModelNode(m_childrenVol[lvIndex], lvName, lvID);
  //   }
  // if (pvCount == 0)
  writePVModelNode(physVol, pvName, pvID);
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

void dd4hep2FBXWriter::writeLVModelNode(Volume logVol, const std::string lvName, unsigned long long lvID)
{
  m_File << "\t; LogVol " << logVol.name() << " with solid " << logVol.solid().name() << std::endl
         << "\tModel: " << lvID << ", \"Model::lv_" << lvName << "\", \"Null\" {" << std::endl
         << "\t\tVersion: 232" << std::endl
         << "\t\tProperties70:  {" << std::endl
         << "\t\t}" << std::endl
         << "\t\tShading: T" << std::endl
         << "\t\tCulling: \"CullingOff\"" << std::endl
         << "\t}" << std::endl;
}

void dd4hep2FBXWriter::writePVModelNode(DetElement physVol, const std::string pvName, unsigned long long pvID)
{
  Position move = physVol.placement().position();
  TGeoRotation mat(physVol.placement().matrix()); //
  // G4RotationMatrix* rot = physVol->GetObjectRotation();
  // G4ThreeVector move = physVol->GetObjectTranslation();
  // FBX uses the Tait-Bryan version of the Euler angles (X then Y then Z rotation)
  double phi = 0, theta = 0, psi = 0;
  mat.GetAngles(phi, theta, psi); // all angles in degrees

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
  m_File << "\t; PhysVol " << physVol.name();
  // if (physVol->IsReplicated()) {
  //   m_File << " (replicated: copy " << physVol->GetCopyNo() << ")";
  // }
  m_File << ", placing LogVol " << physVol.volume().name() << std::endl
         << "\tModel: " << pvID << ", \"Model::" << pvName << "\", \"Null\" {" << std::endl
         << "\t\tVersion: 232" << std::endl
         << "\t\tProperties70:  {" << std::endl
         << "\t\t\tP: \"Lcl Translation\", \"Lcl Translation\", \"\", \"A\"," << move.x() << "," << move.y() << "," << move.z() << std::endl
         << "\t\t\tP: \"Lcl Rotation\", \"Lcl Rotation\", \"\", \"A\"," << psi << "," << theta << "," << psi << std::endl
         << "\t\t}" << std::endl
         << "\t\tShading: T" << std::endl
         << "\t\tCulling: \"CullingOff\"" << std::endl
         << "\t}" << std::endl;
}

void dd4hep2FBXWriter::writeMaterialNode(int lvIndex, const std::string matName)
{
  // G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();
  Volume logVol = m_childrenVol[lvIndex];
  unsigned long long matID = (*m_MatID)[lvIndex];
  float alpha = 0.5, red = 0, green = 1, blue = 0;
  // float color[]= {0.0, 1.0, 0.0,0.5}; // " r, g, b, a// default is semi-transparent green
  // BelleII EM calorimeter crystals
  /*if ((matName.compare(0, 23, "eclBarrelCrystalLogical") == 0) ||
      (matName.compare(0, 20, "eclFwdCrystalLogical") == 0) ||
      (matName.compare(0, 20, "eclBwdCrystalLogical") == 0) ||
      (matName.compare(0, 24, "eclBarrelCrystalPhysical") == 0) ||
      (matName.compare(0, 21, "eclFwdCrystalPhysical") == 0) ||
      (matName.compare(0, 21, "eclBwdCrystalPhysical") == 0)) {
    color = G4Color(1.0, 1.0, 1.0, 1.0); // white since ECL crystals have no G4VisAttribute :(
  }*/
  bool visible = true;
  string materialName = logVol.material().name();
  // Hide volumes that contain vacuum, air or gas
  // if (materialName == "Vacuum")
  //   visible = false;
  // if (materialName == "G4_AIR")
  //   visible = false;
  // if (materialName == "CDCGas") visible = false; // BelleII
  // if (materialName == "ColdAir") visible = false; // BelleII
  // if (materialName == "STR-DryAir") visible = false; // BelleII
  // if (materialName == "TOPAir") visible = false; // BelleII
  // if (materialName == "TOPVacuum") visible = false; // BelleII
  // Hide volumes that are invisible in the GEANT4 geometry
  const VisAttr visAttr = logVol.visAttributes();
  if (!(visAttr == VisAttr()))
  {
    logVol.visAttributes().argb(alpha, red, green, blue);
    if (!(visAttr.visible()))
      {
        // visible = false;
      }
  }
  else
  {
    // visible = false;
  }
  if (logVol.isSensitive())
    visible = true;
  (*m_Visible)[lvIndex] = visible;
  m_File << "\t; Color for LogVol " << logVol.name() << std::endl
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
  }
  else if (solidtype == "TGeoConeSeg") // ConeSegment
  {
    ConeSegment tmp(solid);
    return new GPolyhedronCons(tmp.rMin1(), tmp.rMax1(), tmp.rMin2(), tmp.rMax2(), tmp.dZ(), tmp.startPhi(), tmp.endPhi());
  }
  else if (solidtype == "TGeoTubeSeg") // Tube___________
  {
    Tube tmp(solid);
    return new GPolyhedronTubs(tmp.rMin(), tmp.rMax(), tmp.dZ(), tmp.startPhi(), tmp.endPhi()) ;
  }
  else if (solidtype == "TGeoCtub") // CutTube
  {
    // https://apc.u-paris.fr/~franco/g4doxy/html/G4CutTubs_8cc-source.html#l01975
  }
  else if (solidtype == "TGeoEltu") // EllipticalTube
  {
    EllipticalTube tmp(solid);
    GPolyhedronTube* eTube = new GPolyhedronTube(0.,1., tmp.dZ());
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
  else if (solidtype == "TGeoArb8")  // EightPointSolid
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