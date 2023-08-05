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
#include "PMT.h"
#include "PMTParam.h"

#include "Math/Vector3D.h"
#include "Math/Transform3D.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;


/// @brief contruct the PMT core Solid. 
/// @param lcdd detector description
/// @param e xml_h
/// @param thicknessIdx the index of <thickness> tag in order which is in xml files. 0 is the first <thickness> tag, 1 is the second ....
/// @return pmt Solid
Solid PMTSolid(Detector &lcdd, xml_h e, int thicknessIdx){
  xml_det_t   x_det    = e;
  string      det_name = x_det.nameStr();
  // int         det_id   = x_det.id();

  vector<string> thincknessName;
  vector<double> thicknessValue;
  for (xml_coll_t xc(x_det, _U(module)); xc; ++xc)
  {
    xml_comp_t xmodule       = xc;
    for(xml_coll_t zc(xmodule, _U(thickness)); zc; ++zc)
    {
      xml_comp_t x_thick = zc;
      thincknessName.push_back(x_thick.nameStr());
      thicknessValue.push_back(x_thick.thickness());
    }
  }

  double thickness  = thicknessValue[thicknessIdx];
  string namePrefix = thincknessName[thicknessIdx];
  PlacedVolume pv;
  Volume pmtVol;
  Solid pmt;

  // iterate the module, contents in module
  for (xml_coll_t xc(x_det, _U(module)); xc; ++xc)
  {
    xml_comp_t xmodule       = xc;

    // construct the ellisolid volume
    xml_comp_t pmtElli = xmodule.child(_U(sphere));
    xml_dim_t  ell_scale    = pmtElli.parameters();
    double     rmin         = pmtElli.rmin();
    double     rmax         = pmtElli.rmax()  + thickness;
    double     xScale       = (ell_scale.x()  + thickness)/rmax;
    double     yScale       = (ell_scale.y()  + thickness)/rmax;
    double     zScale       = (ell_scale.z()  + thickness)/rmax;
 
    Sphere pmtSphere(rmin, rmax);
    Solid pmtElliSolid = Scale(pmtSphere, xScale, yScale, zScale);
    // Volume pmtElliVol("pmtEllisolid", pmtElliSolid, pmtElliMat);
  
    // constuct the polycone volume
    xml_comp_t pmtPolycone  = xmodule.child(_U(polycone));
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

    cout << "@@@ PMTConstructor INFO " << __LINE__ << " One PMT solid name: " << "PMT"+namePrefix << " OK " << endl;

    // pmtVol=Volume("PMTVol", pmt, lcdd.air());
  }

  return pmt;
}


/// @brief 
/// @param lcdd 
/// @param e 
/// @return PMT Volume 
// PlacedVolume PMTVol(Detector &lcdd, xml_h e, SensitiveDetector sens){

//   xml_det_t   x_det    = e;
//   string      det_name = x_det.nameStr();
//   int         det_id   = x_det.id();

//   Solid pmt_solid   = PMTSolid(lcdd, e, 0);
//   Solid body_solid  = PMTSolid(lcdd, e, 1);
//   Solid inner_solid = PMTSolid(lcdd, e, 2);

//   double rEllip=0;
//   double zEllip=0;
//   for (xml_coll_t xc(x_det, _U(module)); xc; ++xc)
//   {
//     xml_comp_t xmodule       = xc; 
//     xml_comp_t pmtElli      = xmodule.child(_U(sphere));
//     xml_dim_t  ell_scale    = pmtElli.parameters();
//     rEllip         = pmtElli.rmax();
//     zEllip         = ell_scale.z() ;
//   }

//   Tube pInnerSep("", 0, rEllip+1E-9*mm, zEllip/2.0 + 1E-9*mm, 0, 2*M_PI);
//   // Volume ttt("", pInnerSep, lcdd.air());

//   Position innerSepDispl(0.,0.,zEllip/2.0 -1E-9*mm);
//   Solid inner1_solid = IntersectionSolid(inner_solid, pInnerSep, innerSepDispl);
//   Solid inner2_solid = SubtractionSolid (inner_solid, pInnerSep, innerSepDispl);

//   Volume pmt_log   ("", pmt_solid,    lcdd.material("Pyrex"));
//   Volume body_log  ("", body_solid,   lcdd.material("Pyrex"));
//   Volume inner1_log("", inner1_solid, lcdd.vacuum()); // half of Ellisolid
//   Volume inner2_log("", inner2_solid, lcdd.vacuum()); // half of Ellisolid + tube

//   PlacedVolume body_phys   = pmt_log.placeVolume(body_log);
//   PlacedVolume inner1_phys = body_log.placeVolume(inner1_log); 
//   PlacedVolume inner2_phys = body_log.placeVolume(inner2_log); 


//   // set the vis attributes
//   inner1_log.setVisAttributes(lcdd.visAttributes("inner1Vis"));
//   inner2_log.setVisAttributes(lcdd.visAttributes("inner2Vis"));
//   body_log.setVisAttributes  (lcdd.visAttributes("bodyVis"));

//   return body_log;

// }

static Ref_t create_detector(Detector &lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t   x_det    = e;
  string      det_name = x_det.nameStr();
  int         det_id   = x_det.id();


  // PlacedVolume body_log=PMTVol(lcdd, e, sens);
  Solid pmt_solid   = PMTSolid(lcdd, e, 0);
  Solid body_solid  = PMTSolid(lcdd, e, 1);
  Solid inner_solid = PMTSolid(lcdd, e, 2);

  double rEllip=0;
  double zEllip=0;
  for (xml_coll_t xc(x_det, _U(module)); xc; ++xc)
  {
    xml_comp_t xmodule       = xc; 
    xml_comp_t pmtElli      = xmodule.child(_U(sphere));
    xml_dim_t  ell_scale    = pmtElli.parameters();
    rEllip         = pmtElli.rmax();
    zEllip         = ell_scale.z() ;
  }

  Tube pInnerSep("", 0, rEllip+1E-9*mm, zEllip/2.0 + 1E-9*mm, 0, 2*M_PI);
  // Volume ttt("", pInnerSep, lcdd.air());

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


  // set the vis attributes
  inner1_log.setVisAttributes(lcdd.visAttributes("inner1Vis"));
  inner2_log.setVisAttributes(lcdd.visAttributes("inner2Vis"));
  body_log.setVisAttributes  (lcdd.visAttributes("bodyVis"));


  DetElement sdet(det_name, det_id);
  Volume motherVol = lcdd.pickMotherVolume(sdet);
  PlacedVolume phv = motherVol.placeVolume(body_log, Position(0, 0, 0));
  phv.addPhysVolID("system", sdet.id()).addPhysVolID("pmt", 1);

  sdet.setPlacement(phv);
  cout << "@@@ PMTConstructor INFO " << __LINE__ << " configure one PMT shape finished " << endl;


  // atouch optical surface to the Volume
  // help:
  // https://geant4-userdoc.web.cern.ch/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html#optical-photon-processes
  // https://geant4-userdoc.web.cern.ch/Doxygen/examples_doc/html/ExampleOpNovice.html
  // https://indico.fnal.gov/event/13730/contributions/21420/attachments/13920/17722/optical.pdf
  // /home/wln/DD4hep_source/examples/OpticalSurfaces/src/OpNovice_geo.cpp
  // https://www.liangye.site/2022/05/12/geant4-optical/
  // Now attach the surface
  // need to fix, wrong parameters
  // OpticalSurfaceManager surfMgr = lcdd.surfaceManager();
  // OpticalSurface waterSurf  = surfMgr.opticalSurface("/world/BubbleDevice#WaterSurface");
  // OpticalSurface airSurf    = surfMgr.opticalSurface("/world/BubbleDevice#AirSurface");

  // BorderSurface  tankSurf   = BorderSurface(lcdd, sdet, "HallTank",   waterSurf, inner2_phys, body_phys);
  // BorderSurface  bubbleSurf = BorderSurface(lcdd, sdet, "TankBubble", airSurf,   inner2_phys, body_phys);
  // bubbleSurf.isValid();
  // tankSurf.isValid();


  // if ( x_slice.isSensitive() )  {
  // if ( 1 )  {
    // sens.setType("tracker");
    // pmt_log.setSensitiveDetector(sens);
    // sensitives.push_back(s_phv);
  // }
  cout << "@@@ PMTConstructor INFO " << __LINE__ << " configure PMT optical surface finished " << endl;





/*

  PMTParam PMT_WP("/home/wln/junosw/data/Detector/Geometry/PMTPos_WP_LPMT.csv", PMTType::Hamamatsu);
  PMTParam PMT_CD_Large("/home/wln/junosw/data/Detector/Geometry/PMTPos_CD_LPMT.csv", "/home/wln/junosw/data/Detector/Geometry/PMTType_CD_LPMT.csv");
  PMTParam PMT_CD_Small("/home/wln/junosw/data/Detector/Geometry/PMTPos_CD_SPMT.csv", "/home/wln/junosw/data/Detector/Geometry/PMTType_CD_SPMT.csv");


  PMTParam test("/home/wln/junosw/data/Detector/Geometry/PMTPos_CD_LPMT.csv", "/home/wln/junosw/data/Detector/Geometry/PMTType_CD_LPMT.csv");
  map<int, PMT> pmtpars = test.getPMTParam();
  for(uint64_t i = 0; i< pmtpars.size();i++)
  {
    if(i>10) break;
    cout << pmtpars[i].getPMTPosition().X() << " " << 
    pmtpars[i].getPMTPosition().Y() << " " << 
    pmtpars[i].getPMTPosition().X() << " " << " "
     << int(pmtpars[i].getPMTType()) << " " << endl;

  }
  
  cout << "@@ " <<pmtpars.size() <<endl;
  cout << 
  int(PMTType::NotSet) << " " <<
  int(PMTType::HZC) <<" " << 
  int(PMTType::Hamamatsu) <<" " << 
  int(PMTType::NNVT) << " " <<
  int(PMTType::HighQENNVT) << " " <<
  
  endl;


  // {NotSet=1, HZC, Hamamatsu, NNVT, HighQENNVT};

  */

  







  return sdet;
}

DECLARE_DETELEMENT(DD4hep_CentralDetector_PMT, create_detector)
