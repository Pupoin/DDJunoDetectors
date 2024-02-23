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

#include "GPolyhedron.hh"

GPolyhedron::GPolyhedron ():
  fNumberOfRotationStepsAtTimeOfCreation (fNumberOfRotationSteps)
{}

GPolyhedron::~GPolyhedron () {}

GPolyhedron::GPolyhedron (const HepPolyhedron& from)
  : HepPolyhedron(from)
{
  fNumberOfRotationStepsAtTimeOfCreation =
    from.fNumberOfRotationSteps;
}

GPolyhedronBox::GPolyhedronBox (double dx, double dy, double dz):
  GPolyhedron (HepPolyhedronBox (dx, dy, dz)) {}

GPolyhedronBox::~GPolyhedronBox () {}

GPolyhedronCone::GPolyhedronCone (double Rmn1, double Rmx1,
                                    double Rmn2, double Rmx2, double Dz):
  GPolyhedron (HepPolyhedronCone (Rmn1, Rmx1, Rmn2, Rmx2, Dz)) {}

GPolyhedronCone::~GPolyhedronCone () {}

GPolyhedronCons::GPolyhedronCons (double Rmn1, double Rmx1,
                                    double Rmn2, double Rmx2, double Dz,
                                    double Phi1, double Dphi):
  GPolyhedron (HepPolyhedronCons (Rmn1, Rmx1, Rmn2, Rmx2, Dz, Phi1, Dphi)) {}

GPolyhedronCons::~GPolyhedronCons () {}

GPolyhedronPara::GPolyhedronPara (double Dx, double Dy, double Dz,
                                    double Alpha, double Theta,
                                    double Phi):
  GPolyhedron (HepPolyhedronPara (Dx, Dy, Dz, Alpha, Theta, Phi)) {}

GPolyhedronPara::~GPolyhedronPara () {}

GPolyhedronPcon::GPolyhedronPcon (double phi, double dphi, int nz,
                                    const double *z,
                                    const double *rmin,
                                    const double *rmax):
  GPolyhedron (HepPolyhedronPcon (phi, dphi, nz, z, rmin, rmax)) {}

GPolyhedronPcon::GPolyhedronPcon (double phi, double dphi,
                                    const std::vector<CLHEP::Hep2Vector> &rz):
  GPolyhedron (HepPolyhedronPcon(phi, dphi, rz)) {}

GPolyhedronPcon::~GPolyhedronPcon () {}

GPolyhedronPgon::GPolyhedronPgon (double phi, double dphi, int npdv,
                                    int nz,
                                    const double *z,
                                    const double *rmin,
                                    const double *rmax):
  GPolyhedron (HepPolyhedronPgon (phi, dphi, npdv, nz, z, rmin, rmax)) {}

GPolyhedronPgon::GPolyhedronPgon (double phi, double dphi, int npdv,
                                    const std::vector<CLHEP::Hep2Vector> &rz):
  GPolyhedron (HepPolyhedronPgon(phi, dphi, npdv, rz)) {}

GPolyhedronPgon::~GPolyhedronPgon () {}

GPolyhedronSphere::GPolyhedronSphere (double rmin, double rmax,
                                        double phi, double dphi,
                                        double the, double dthe):
  GPolyhedron (HepPolyhedronSphere (rmin, rmax, phi, dphi, the, dthe)) {}

GPolyhedronSphere::~GPolyhedronSphere () {}

GPolyhedronTet::GPolyhedronTet (const double p0[3],
                                  const double p1[3],
                                  const double p2[3],
                                  const double p3[3]):
  GPolyhedron (HepPolyhedronTet (p0, p1, p2, p3)) {}

GPolyhedronTet::~GPolyhedronTet () {}

GPolyhedronTorus::GPolyhedronTorus (double rmin, double rmax,
                                      double rtor,
                                      double phi, double dphi):
  GPolyhedron (HepPolyhedronTorus (rmin, rmax, rtor, phi, dphi)) {}

GPolyhedronTorus::~GPolyhedronTorus () {}

GPolyhedronTrap::GPolyhedronTrap (double Dz, double Theta, double Phi,
                                    double Dy1,
                                    double Dx1, double Dx2, double Alp1,
                                    double Dy2,
                                    double Dx3, double Dx4, double Alp2):
  GPolyhedron (HepPolyhedronTrap (Dz, Theta, Phi, Dy1, Dx1, Dx2, Alp1,
                                   Dy2, Dx3, Dx4, Alp2)) {}

GPolyhedronTrap::~GPolyhedronTrap () {}

GPolyhedronTrd1::GPolyhedronTrd1 (double Dx1, double Dx2,
                                    double Dy, double Dz):
  GPolyhedron (HepPolyhedronTrd1 (Dx1, Dx2, Dy, Dz)) {}

GPolyhedronTrd1::~GPolyhedronTrd1 () {}

GPolyhedronTrd2::GPolyhedronTrd2 (double Dx1, double Dx2,
                                    double Dy1, double Dy2, double Dz):
  GPolyhedron (HepPolyhedronTrd2 (Dx1, Dx2, Dy1, Dy2, Dz)) {}

GPolyhedronTrd2::~GPolyhedronTrd2 () {}

GPolyhedronTube::GPolyhedronTube (double Rmin, double Rmax, double Dz):
  GPolyhedron (HepPolyhedronTube (Rmin, Rmax, Dz)) {}

GPolyhedronTube::~GPolyhedronTube () {}

GPolyhedronTubs::GPolyhedronTubs (double Rmin, double Rmax, double Dz,
                                    double Phi1, double Dphi):
  GPolyhedron (HepPolyhedronTubs (Rmin, Rmax, Dz, Phi1, Dphi)) {}

GPolyhedronTubs::~GPolyhedronTubs () {}

GPolyhedronParaboloid::GPolyhedronParaboloid (double r1, double r2,
                                                double dz, double sPhi,
                                                double dPhi):
  GPolyhedron (HepPolyhedronParaboloid(r1, r2, dz, sPhi, dPhi)) {}

GPolyhedronParaboloid::~GPolyhedronParaboloid () {}

GPolyhedronHype::GPolyhedronHype (double r1, double r2, double tan1,
                                    double tan2, double halfZ):
  GPolyhedron (HepPolyhedronHype(r1, r2, tan1, tan2, halfZ)) {}

GPolyhedronHype::~GPolyhedronHype () {}

GPolyhedronEllipsoid::GPolyhedronEllipsoid (double ax, double by,
                                              double cz,
                                              double zCut1, double zCut2):
  GPolyhedron (HepPolyhedronEllipsoid (ax, by, cz, zCut1, zCut2)) {}

GPolyhedronEllipsoid::~GPolyhedronEllipsoid () {}

GPolyhedronEllipticalCone::GPolyhedronEllipticalCone (double ax,
                                                        double ay,
                                                        double h,
                                                        double zCut1):
  GPolyhedron (HepPolyhedronEllipticalCone (ax, ay, h, zCut1)) {}

GPolyhedronEllipticalCone::~GPolyhedronEllipticalCone () {}

GPolyhedronHyperbolicMirror::GPolyhedronHyperbolicMirror (double a,
                                                            double h,
                                                            double r):
  GPolyhedron (HepPolyhedronHyperbolicMirror(a, h, r)) {}

GPolyhedronHyperbolicMirror::~GPolyhedronHyperbolicMirror () {}

std::ostream& operator<<(std::ostream& os, const GPolyhedron& polyhedron)
{
  os << "GPolyhedron: "
    //  << (const G4Visible&)polyhedron << '\n'
     << (const HepPolyhedron&)polyhedron;
  return os;
}
