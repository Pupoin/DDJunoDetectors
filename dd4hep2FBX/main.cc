#include "FBXWriter.h"

#include "G4GDMLParser.hh"
#include "G4VPhysicalVolume.hh"

#include <iomanip>
#include <iostream>
#include <string>

int main()
{
  std::cout << "Hello World!" << std::endl;

  G4GDMLParser parser;
  //parser.Read("MucTest.gdml");
  //parser.Read("muonc.gdml");
  //parser.Read("Emc.gdml");
  parser.Read("EicC.gdml");


  G4VPhysicalVolume* fWorld = 0;
  fWorld = (G4VPhysicalVolume *)parser.GetWorldVolume();
  if (0 == fWorld) {
    std::cout << "World volume not found!" << std::endl;
  }

  FBXWriter b;
  b.doit("out.fbx", true);
  std::cout << "End of FBXWriter" << std::endl;
}