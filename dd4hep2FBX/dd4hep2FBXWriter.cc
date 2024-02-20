/**************************************************************************
 * BASF2 (Belle Analysis Framework 2)                                     *
 * Copyright(C) 2016 - Belle II Collaboration                             *
 *                                                                        *
 * Author: ZhaoYang Yuan                                                   *
 * 2024/02/20                                                                        *
 * This software is provided "as is" without any warranty.                *
 **************************************************************************/

// #include <dd4hep2FBXWriter.h>

#include "DD4hep/DetFactoryHelper.h"
#include "XML/Layering.h"

#include "Math/Vector3D.h"
#include "Math/Transform3D.h"

using namespace std;
using namespace dd4hep;
using namespace dd4hep::detail;

#include <iomanip>

int main(){
  Detector & lcwdd = Detector::getInstance();
  Tube waterPool(0, 10*mm, 10*mm, 20*mm, 30*mm);
  lcwdd.fromCompact("/home/wln/DD4hep/DDJunoDetectors/compact/Juno.xml");
  // lcwdd.fromCompact("/home/wln/DD4hep2FBX/lcdd/watePool.xml");
  DetElement ele=lcwdd.detector("WaterPool");
  std::cout << "@@@@@ " << ele.id() << "  " << ele.name() << std::endl;

  return 0;
}