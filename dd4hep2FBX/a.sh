#!/bin/bash
 
#include <CLHEP/Vector/TwoVector.h> //#include "G4TwoVector.hh"
#include <CLHEP/Geometry/Point3D.h>  //#include "G4Point3D.hh"
#include <CLHEP/Geometry/Normal3D.h> //#include "G4Normal3D.hh"
#include <CLHEP/Geometry/Transform3D.h> //#include "G4Transform3D.hh"
#include <CLHEP/Geometry/Vector3D.h> //#include "G4Vector3D.hh"

ips=(
./src/BooleanProcessor.src
./src/G4Polyhedron.cc
./src/HepPolyhedron.cc
./src/HepPolyhedronProcessor.src
./include/G4Polyhedron.hh
./include/HepPolyhedron.h
./include/HepPolyhedronProcessor.h
# ./a
)
for ip in ${ips[@]}
do
    sed "s|G4int|int|g" -i ${ip}
    sed "s|G4double|double|g" -i ${ip}
    sed "s|G4cout|std::cout|g" -i ${ip}
    sed "s|G4endl|std::endl|g" -i ${ip}
    sed "s|G4cerr|std::cerr|g" -i ${ip}
    sed "s|G4bool|bool|g" -i ${ip}

    # sed "s|G4TwoVector|CLHEP::Hep2Vector|g" -i ${ip}           
    # sed "s|G4Point3D|HepGeom::Point3D<double>|g" -i ${ip} 
    # sed "s|G4Normal3D|HepGeom::Normal3D<double>|g" -i ${ip} 
    # sed "s|G4Transform3D|HepGeom::Transform3D|g" -i ${ip} 	
    # break;
done
