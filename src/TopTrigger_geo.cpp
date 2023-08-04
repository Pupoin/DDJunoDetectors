//==========================================================================
//  AIDA Detector description implementation
//--------------------------------------------------------------------------
// Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
// All rights reserved.
//
// For the licensing terms see $DD4hepINSTALL/LICENSE.
// For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
//
// Author     : Zhaoyang Yuan
// Data       : 2023-08-01
//
//==========================================================================
//
// Specialized Top Trigger detector constructor
//
//==========================================================================
#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;
static Ref_t create_detector(Detector &lcdd, xml_h e, SensitiveDetector sens)
{
  xml_det_t   x_det    = e;
  string      det_name = x_det.nameStr();
  int         det_id   = x_det.id();

  PlacedVolume pv;

  Volume TopTriggerVol("", Box(1,1,1), lcdd.air());

  // iterate the layer, contents(spheres) in layer
  for (xml_coll_t xlayer(x_det, _U(layer)); xlayer; ++xlayer)
  {
    for (xml_coll_t xc(xlayer, _U(sphere)); xc; ++xc)
    {
      xml_comp_t mmbox     = xc;
      string     s_name    = mmbox.nameStr();
      double    box_x      = mmbox.x();
      double    box_y      = mmbox.y();
      double    box_z      = mmbox.z();

      xml_comp_t mmposition  = mmbox.position();
      double    pos_x        = mmposition.x();
      double    pos_y        = mmposition.y();
      double    pos_z        = mmposition.z();

      Box       tt(box_x, box_y, box_z);
      Volume ttVol("", tt, lcdd.air());
      TopTriggerVol.placeVolume(ttVol,Position(pos_x,pos_y,pos_z));

      if (mmbox.isSensitive())
      {
        sens.setType(s_name);
        ttVol.setSensitiveDetector(sens);
        // sensitives.push_back(s_phv);
      }
      cout << "TopTrigger_geo INFO @@@ " << __LINE__ << " set " << s_name << " " << "" << endl;
    }
  }

  
  DetElement sdet(det_name, det_id);
  Volume motherVol = lcdd.pickMotherVolume(sdet);
  PlacedVolume phv = motherVol.placeVolume(TopTriggerVol, Position(0, 0, 0));
  // phv.addPhysVolID("system", sdet.id()).addPhysVolID("sphere", 1);

  sdet.setPlacement(phv);
  cout << "TopTrigger_geo INFO @@@ " << __LINE__ << " configure total Top trigger finished " << endl;

  return sdet;
}

DECLARE_DETELEMENT(DD4hep_TopTrigger, create_detector)
