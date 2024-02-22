//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

#include "G4Polyhedron.hh"

G4Polyhedron::G4Polyhedron ():
  fNumberOfRotationStepsAtTimeOfCreation (fNumberOfRotationSteps)
{}

G4Polyhedron::~G4Polyhedron () {}

G4Polyhedron::G4Polyhedron (const HepPolyhedron& from)
  : HepPolyhedron(from)
{
  fNumberOfRotationStepsAtTimeOfCreation =
    from.fNumberOfRotationSteps;
}

G4PolyhedronBox::G4PolyhedronBox (double dx, double dy, double dz):
  G4Polyhedron (HepPolyhedronBox (dx, dy, dz)) {}

G4PolyhedronBox::~G4PolyhedronBox () {}

G4PolyhedronCone::G4PolyhedronCone (double Rmn1, double Rmx1,
                                    double Rmn2, double Rmx2, double Dz):
  G4Polyhedron (HepPolyhedronCone (Rmn1, Rmx1, Rmn2, Rmx2, Dz)) {}

G4PolyhedronCone::~G4PolyhedronCone () {}

G4PolyhedronCons::G4PolyhedronCons (double Rmn1, double Rmx1,
                                    double Rmn2, double Rmx2, double Dz,
                                    double Phi1, double Dphi):
  G4Polyhedron (HepPolyhedronCons (Rmn1, Rmx1, Rmn2, Rmx2, Dz, Phi1, Dphi)) {}

G4PolyhedronCons::~G4PolyhedronCons () {}

G4PolyhedronPara::G4PolyhedronPara (double Dx, double Dy, double Dz,
                                    double Alpha, double Theta,
                                    double Phi):
  G4Polyhedron (HepPolyhedronPara (Dx, Dy, Dz, Alpha, Theta, Phi)) {}

G4PolyhedronPara::~G4PolyhedronPara () {}

G4PolyhedronPcon::G4PolyhedronPcon (double phi, double dphi, int nz,
                                    const double *z,
                                    const double *rmin,
                                    const double *rmax):
  G4Polyhedron (HepPolyhedronPcon (phi, dphi, nz, z, rmin, rmax)) {}

G4PolyhedronPcon::G4PolyhedronPcon (double phi, double dphi,
                                    const std::vector<CLHEP::Hep2Vector> &rz):
  G4Polyhedron (HepPolyhedronPcon(phi, dphi, rz)) {}

G4PolyhedronPcon::~G4PolyhedronPcon () {}

G4PolyhedronPgon::G4PolyhedronPgon (double phi, double dphi, int npdv,
                                    int nz,
                                    const double *z,
                                    const double *rmin,
                                    const double *rmax):
  G4Polyhedron (HepPolyhedronPgon (phi, dphi, npdv, nz, z, rmin, rmax)) {}

G4PolyhedronPgon::G4PolyhedronPgon (double phi, double dphi, int npdv,
                                    const std::vector<CLHEP::Hep2Vector> &rz):
  G4Polyhedron (HepPolyhedronPgon(phi, dphi, npdv, rz)) {}

G4PolyhedronPgon::~G4PolyhedronPgon () {}

G4PolyhedronSphere::G4PolyhedronSphere (double rmin, double rmax,
                                        double phi, double dphi,
                                        double the, double dthe):
  G4Polyhedron (HepPolyhedronSphere (rmin, rmax, phi, dphi, the, dthe)) {}

G4PolyhedronSphere::~G4PolyhedronSphere () {}

G4PolyhedronTet::G4PolyhedronTet (const double p0[3],
                                  const double p1[3],
                                  const double p2[3],
                                  const double p3[3]):
  G4Polyhedron (HepPolyhedronTet (p0, p1, p2, p3)) {}

G4PolyhedronTet::~G4PolyhedronTet () {}

G4PolyhedronTorus::G4PolyhedronTorus (double rmin, double rmax,
                                      double rtor,
                                      double phi, double dphi):
  G4Polyhedron (HepPolyhedronTorus (rmin, rmax, rtor, phi, dphi)) {}

G4PolyhedronTorus::~G4PolyhedronTorus () {}

G4PolyhedronTrap::G4PolyhedronTrap (double Dz, double Theta, double Phi,
                                    double Dy1,
                                    double Dx1, double Dx2, double Alp1,
                                    double Dy2,
                                    double Dx3, double Dx4, double Alp2):
  G4Polyhedron (HepPolyhedronTrap (Dz, Theta, Phi, Dy1, Dx1, Dx2, Alp1,
                                   Dy2, Dx3, Dx4, Alp2)) {}

G4PolyhedronTrap::~G4PolyhedronTrap () {}

G4PolyhedronTrd1::G4PolyhedronTrd1 (double Dx1, double Dx2,
                                    double Dy, double Dz):
  G4Polyhedron (HepPolyhedronTrd1 (Dx1, Dx2, Dy, Dz)) {}

G4PolyhedronTrd1::~G4PolyhedronTrd1 () {}

G4PolyhedronTrd2::G4PolyhedronTrd2 (double Dx1, double Dx2,
                                    double Dy1, double Dy2, double Dz):
  G4Polyhedron (HepPolyhedronTrd2 (Dx1, Dx2, Dy1, Dy2, Dz)) {}

G4PolyhedronTrd2::~G4PolyhedronTrd2 () {}

G4PolyhedronTube::G4PolyhedronTube (double Rmin, double Rmax, double Dz):
  G4Polyhedron (HepPolyhedronTube (Rmin, Rmax, Dz)) {}

G4PolyhedronTube::~G4PolyhedronTube () {}

G4PolyhedronTubs::G4PolyhedronTubs (double Rmin, double Rmax, double Dz,
                                    double Phi1, double Dphi):
  G4Polyhedron (HepPolyhedronTubs (Rmin, Rmax, Dz, Phi1, Dphi)) {}

G4PolyhedronTubs::~G4PolyhedronTubs () {}

G4PolyhedronParaboloid::G4PolyhedronParaboloid (double r1, double r2,
                                                double dz, double sPhi,
                                                double dPhi):
  G4Polyhedron (HepPolyhedronParaboloid(r1, r2, dz, sPhi, dPhi)) {}

G4PolyhedronParaboloid::~G4PolyhedronParaboloid () {}

G4PolyhedronHype::G4PolyhedronHype (double r1, double r2, double tan1,
                                    double tan2, double halfZ):
  G4Polyhedron (HepPolyhedronHype(r1, r2, tan1, tan2, halfZ)) {}

G4PolyhedronHype::~G4PolyhedronHype () {}

G4PolyhedronEllipsoid::G4PolyhedronEllipsoid (double ax, double by,
                                              double cz,
                                              double zCut1, double zCut2):
  G4Polyhedron (HepPolyhedronEllipsoid (ax, by, cz, zCut1, zCut2)) {}

G4PolyhedronEllipsoid::~G4PolyhedronEllipsoid () {}

G4PolyhedronEllipticalCone::G4PolyhedronEllipticalCone (double ax,
                                                        double ay,
                                                        double h,
                                                        double zCut1):
  G4Polyhedron (HepPolyhedronEllipticalCone (ax, ay, h, zCut1)) {}

G4PolyhedronEllipticalCone::~G4PolyhedronEllipticalCone () {}

G4PolyhedronHyperbolicMirror::G4PolyhedronHyperbolicMirror (double a,
                                                            double h,
                                                            double r):
  G4Polyhedron (HepPolyhedronHyperbolicMirror(a, h, r)) {}

G4PolyhedronHyperbolicMirror::~G4PolyhedronHyperbolicMirror () {}

std::ostream& operator<<(std::ostream& os, const G4Polyhedron& polyhedron)
{
  os << "G4Polyhedron: "
    //  << (const G4Visible&)polyhedron << '\n'
     << (const HepPolyhedron&)polyhedron;
  return os;
}
