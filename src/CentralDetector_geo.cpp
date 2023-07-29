//==========================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : M.Frank
//
//==========================================================================
//
// Specialized generic detector constructor
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
static Ref_t create_detector(Detector &lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t x_det = e;
  string det_name = x_det.nameStr();
  int det_id = x_det.id();

  DetElement sdet(det_name, det_id);
  Volume motherVol = lcdd.pickMotherVolume(sdet);
  PlacedVolume pv;



  /// get the minimum(rmin) and maximum(rmax) of all spheres
  double totRmin = 999999 * dd4hep::m, totRmax = 0 * dd4hep::m;
  for (xml_coll_t xlayer(x_det, _U(layer)); xlayer; ++xlayer)
  {
    for (xml_coll_t xc(xlayer, _U(sphere)); xc; ++xc)
    {
      xml_comp_t x_sphere = xc;
      double rmin = x_sphere.rmin();
      double rmax = x_sphere.rmax();

      if (rmin <= totRmin)
        totRmin = rmin;
      if (rmax >= totRmax)
        totRmax = rmax;
    }
  }
  
  Material air = lcdd.air();
  Sphere totSphere(totRmin, totRmax, 0, M_PI, 0, 2 * M_PI);
  Volume totSphereVol(det_name + "_totSphere", totSphere, air);

  // iterate the layer, contents(spheres) in layer
  for (xml_coll_t xlayer(x_det, _U(layer)); xlayer; ++xlayer)
  {
    int sphereNum = 0;
    for (xml_coll_t xc(xlayer, _U(sphere)); xc; ++xc)
    {
      xml_comp_t x_sphere = xc;
      double rmin = x_sphere.rmin();
      double rmax = x_sphere.rmax();
      double startTheta = x_sphere.starttheta();
      double endTheta = x_sphere.endtheta();
      double startPhi = x_sphere.startphi();
      double endPhi = x_sphere.endphi();
      Material materi = lcdd.material(x_sphere.materialStr());
      string s_name = x_sphere.nameStr();
      double pos_x = x_sphere.x();
      double pos_y = x_sphere.y();
      double pos_z = x_sphere.z();

      sphereNum++;

      Sphere LiquidScintillator(rmin, rmax, startTheta, endTheta, startPhi, endPhi);
      Volume LSVol(s_name, LiquidScintillator, materi);
      LSVol.setVisAttributes(lcdd.visAttributes( x_sphere.visStr()));
      pv = totSphereVol.placeVolume(LSVol, Position(pos_x, pos_y, pos_z));
      // pv = motherVol.placeVolume(v_LS, Position(0, 0, 0));
      pv.addPhysVolID("sphere", sphereNum);

      if (x_sphere.isSensitive())
      {
        sens.setType(s_name);
        LSVol.setSensitiveDetector(sens);
        // sensitives.push_back(s_phv);
      }
      cout << "TotSphere_geo INFO @@@ " << __LINE__ << " set " << s_name << " " << sphereNum << endl;
    }
  }

  

  PlacedVolume phv = motherVol.placeVolume(totSphereVol, Position(0, 0, 0));
  phv.addPhysVolID("system", sdet.id()).addPhysVolID("sphere", 1);

  sdet.setPlacement(phv);
  cout << "TotSphere_geo INFO @@@ " << __LINE__ << " configure total spheres finished " << endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_CentralDetector, create_detector)
