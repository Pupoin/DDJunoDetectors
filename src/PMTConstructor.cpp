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


/// @brief contruct the PMT core Solid. 
/// @param lcdd detector description
/// @param e xml_h
/// @param sens SensitiveDetector
/// @param thicknessIdx the index of <thickness> tag in order which is in xml files. 0 is the first <thickness> tag, 1 is the second ....
/// @return Solid
Solid PMTCore(Detector &lcdd, xml_h e, int thicknessIdx){
  xml_det_t   x_det    = e;
  string      det_name = x_det.nameStr();
  // int         det_id   = x_det.id();

  vector<string> thincknessName;
  vector<double> thicknessValue;
  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc)
  {
    xml_comp_t xlayer       = xc;
    for(xml_coll_t zc(xlayer, _U(thickness)); zc; ++zc)
    {
      xml_comp_t x_thick = zc;
      thincknessName.push_back(x_thick.nameStr());
      thicknessValue.push_back(x_thick.thickness());
    }
  }

  double thickness  = thicknessValue[thicknessIdx];
  string namePrefix = thincknessName[thicknessIdx];
  cout <<"fafasdadffff " <<thickness << endl;
  PlacedVolume pv;
  Volume pmtVol;
  Solid pmt;

  // iterate the layer, contents in layer
  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc)
  {
    xml_comp_t xlayer       = xc;

    // construct the ellisolid volume
    xml_comp_t pmtElli = xlayer.child(_U(sphere));
    // xml_dim_t  ell_position = pmtElli.position();
    xml_dim_t  ell_scale    = pmtElli.parameters();
    // Material   pmtElliMat    = lcdd.material(pmtElli.materialStr());
    double     rmin         = pmtElli.rmin();
    double     rmax         = pmtElli.rmax()  + thickness;
    double     xScale       = (ell_scale.x()  + thickness)/rmax;
    double     yScale       = (ell_scale.y()  + thickness)/rmax;
    double     zScale       = (ell_scale.z()  + thickness)/rmax;
 
    Sphere pmtSphere(rmin, rmax);
    Solid pmtElliSolid = Scale(pmtSphere, xScale, yScale, zScale);
    // Volume pmtElliVol("pmtEllisolid", pmtElliSolid, pmtElliMat);
  
    // constuct the polycone volume
    xml_comp_t pmtPolycone  = xlayer.child(_U(polycone));
    // Material   pmtPolMat    = lcdd.material(pmtPolycone.materialStr()); 
    double     pol_startPhi = pmtPolycone.start();
    double     pol_deltaPhi = pmtPolycone.deltaphi();
    vector<double> zPlane, rInner, rOuter;
    vector<double> zplane_z={-thickness, thickness, thickness, 0};
    int idx=0;
    for(xml_coll_t zc(pmtPolycone, _U(zplane)); zc; ++zc)
    {
        xml_comp_t x_zplane = zc;
        zPlane.push_back(x_zplane.z() + zplane_z[idx]);
        rInner.push_back(x_zplane.rmin());
        rOuter.push_back(x_zplane.rmax()+thickness);
        idx++;
    }
    Solid pmtPolSolid = Polycone(pol_startPhi, pol_deltaPhi, rInner, rOuter, zPlane);
    // Volume pmtPolycVol("pmtPolycone", pmtPolSolid, pmtPolMat);


    // ellisolid + polycone
    Transform3D tr(EulerAngles(0,0,0), Translation3D (0,0,0));
    pmt = UnionSolid("PMT"+namePrefix, pmtElliSolid, pmtPolSolid, tr);

    cout << "PMTConstructor INFO @@@ " << __LINE__ << " One PMT solid name: " << "PMT"+namePrefix << " OK " << endl;

    // pmtVol=Volume("PMTVol", pmt, lcdd.air());
  }

  return pmt;
}


static Ref_t create_detector(Detector &lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t   x_det    = e;
  string      det_name = x_det.nameStr();
  int         det_id   = x_det.id();

  Solid pmt_solid   = PMTCore(lcdd, e, 0);
  Solid body_solid  = PMTCore(lcdd, e, 1);
  Solid inner_solid = PMTCore(lcdd, e, 2);

  double rEllip=0;
  double zEllip=0;
  for (xml_coll_t xc(x_det, _U(layer)); xc; ++xc)
  {
    xml_comp_t xlayer       = xc;
    xml_comp_t pmtElli      = xlayer.child(_U(sphere));
    xml_dim_t  ell_scale    = pmtElli.parameters();
    rEllip         = pmtElli.rmax();
    zEllip         = ell_scale.z() ;
  }

  Tube pInnerSep("", 0, rEllip+1E-9*mm, zEllip/2.0 + 1E-9*mm, 0, 2*M_PI);
  Volume ttt("", pInnerSep, lcdd.air());

  Position innerSepDispl(0.,0.,zEllip/2.0 -1E-9*mm);
  Solid inner1_solid = IntersectionSolid(inner_solid, pInnerSep, innerSepDispl);
  Solid inner2_solid = SubtractionSolid (inner_solid, pInnerSep, innerSepDispl);

  Volume pmt_log   ("", pmt_solid,    lcdd.material("Pyrex"));
  Volume body_log  ("", body_solid,   lcdd.material("Pyrex"));
  Volume inner1_log("", inner1_solid, lcdd.vacuum()); // half of Ellisolid
  Volume inner2_log("", inner2_solid, lcdd.vacuum()); // half of Ellisolid + tube

  PlacedVolume body_phys   = pmt_log.placeVolume(body_log);
  PlacedVolume inner1_phys = body_log.placeVolume(inner1_log); 
  PlacedVolume inner2_phys = body_log.placeVolume(inner2_log); 


  inner1_log.setVisAttributes(lcdd.visAttributes("inner1Vis"));
  inner2_log.setVisAttributes(lcdd.visAttributes("inner2Vis"));
  body_log.setVisAttributes(lcdd.visAttributes("bodyVis"));




  DetElement sdet(det_name, det_id);
  Volume motherVol = lcdd.pickMotherVolume(sdet);
  PlacedVolume phv = motherVol.placeVolume(body_log, Position(0, 0, 0));

  // phv.addPhysVolID("system", sdet.id()).addPhysVolID("sphere", 1);

  sdet.setPlacement(phv);
  cout << "TotSphere_geo INFO @@@ " << __LINE__ << " configure total spheres finished " << endl;

  




  return sdet;
}

DECLARE_DETELEMENT(DD4hep_CentralDetector_PMT, create_detector)