//==========================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : Zhaoyang Yuan 2023-07-30
//
//==========================================================================
//
// Specialized generic detector constructor
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"

#include "Math/Vector3D.h"
#include "Math/Transform3D.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

static Ref_t create_detector(Detector &lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det    = e;
  string    det_name = x_det.nameStr();
  int       det_id   = x_det.id();
  xml_comp_t  mmtube = x_det.child(_U(tube));
  double rmin     = mmtube.rmin();
  double rmax     = mmtube.rmax();
  double startphi = mmtube.startphi();
  double endphi   = mmtube.endphi();
  double height   = mmtube.height();

  Tube waterPool(rmin, rmax, height/2, startphi, endphi);
  Volume wpVol("", waterPool, lcdd.material(mmtube.materialStr()));
  // G_WATERPOOLVOL = Volume("WaterPool", waterPool, lcdd.material(mmtube.materialStr()));
  wpVol.setVisAttributes(lcdd.visAttributes(mmtube.visStr()));

  DetElement sdet(det_name, det_id);
  Volume motherVol = lcdd.pickMotherVolume(sdet);
  PlacedVolume phv = motherVol.placeVolume(wpVol, Position(0, 0, 0));

  // phv.addPhysVolID("system", sdet.id()).addPhysVolID("sphere", 1);

  sdet.setPlacement(phv);
  cout << "WaterPool_geo INFO @@@ " << __LINE__ << " configure water pool finished " << endl;

  return sdet;
  // return NULL;
}

DECLARE_DETELEMENT(DD4hep_WaterPool, create_detector)
