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

#ifndef GPolyhedron_HH
#define GPolyhedron_HH

// Class Description:
// GPolyhedron is an intermediate class between G4 and visualization
// systems. It is intended to provide some service like:
//   - polygonization of the G4 shapes with triangulization
//     (quadrilaterization) of complex polygons;
//   - calculation of normals for faces and vertices.
//
// Inherits from HepPolyhedron, to which reference should be made for
// functionality.
//
// Public constructors:
//   GPolyhedronBox(dx,dy,dz)            - create GPolyhedron for G4 Box;
//   GPolyhedronTrd1(dx1,dx2,dy,dz)      - create GPolyhedron for G4 Trd1;
//   GPolyhedronTrd2(dx1,dx2,dy1,dy2,dz) - create GPolyhedron for G4 Trd2;
//   GPolyhedronTrap(dz,theta,phi,
//                    h1,bl1,tl1,alp1,
//                    h2,bl2,tl2,alp2)    - create GPolyhedron for G4 Trap;
//   GPolyhedronPara(dx,dy,dz,
//                    alpha,theta,phi)    - create GPolyhedron for G4 Para;
//
//   GPolyhedronTube(rmin,rmax,dz)       - create GPolyhedron for G4 Tube;
//   GPolyhedronTubs(rmin,rmax,dz,
//                    phi1,dphi)          - create GPolyhedron for G4 Tubs;
//   GPolyhedronCone(rmin1,rmax1,
//                    rmin2,rmax2,dz)     - create GPolyhedron for G4 Cone;
//   GPolyhedronCons(rmin1,rmax1,
//                    rmin2,rmax2,dz,
//                    phi1,dphi)          - create GPolyhedron for G4 Cons;
//
//   GPolyhedronPgon(phi,dphi,npdv,nz,
//                    z(*),rmin(*),rmax(*)) - create GPolyhedron for G4 Pgon;
//   GPolyhedronPcon(phi,dphi,nz,
//                    z(*),rmin(*),rmax(*)) - create GPolyhedron for G4 Pcon;
//
//   GPolyhedronSphere(rmin,rmax,
//                      phi,dphi,the,dthe)  - create GPolyhedron for Sphere;
//   GPolyhedronTorus(rmin,rmax,rtor,
//                     phi,dphi)            - create GPolyhedron for Torus;
//   GPolyhedronTet(p0[3],p1[3],p2[3],p3[3]) - create polyhedron for Tet;
//
//   GPolyhedronEllipsoid(dx,dy,dz,
//                         zcut1,zcut2)     - create GPolyhedron for Ellipsoid;
//   GPolyhedronEllipticalCone(dx,dy,z,
//                              zcut1)      - create polyhedron for Elliptical cone;
//   GPolyhedronParaboloid(r1,r2,dz,
//                          phi,dphi)       - create polyhedron for Paraboloid;
//   GPolyhedronHype(r1,r2,
//                    tan1,tan2,halfz)      - create polyhedron for Hype;
//   GPolyhedronHyperbolicMirror(a,h,r)    - create polyhedron for Hyperbolic mirror;
//
// Public functions inherited from HepPolyhedron (this list might be
// incomplete):
//   GetNoVertices()  - returns number of vertices
//   GetNoFacets()    - returns number of faces
//   GetNextVertexIndex(index, edgeFlag) - get vertex indeces of the
//                      quadrilaterals in order; returns false when
//                      finished each face;
//   GetVertex(index) - returns vertex by index;
//   GetNextVertex(vertex, edgeFlag) - get vertices with edge visibility
//                      of the quadrilaterals in order;
//                      returns false when finished each face;
//   GetNextVertex(vertex, edgeFlag, normal) - get vertices with edge
//                      visibility and normal of the quadrilaterals
//                      in order; returns false when finished each face;
//   GetNextNormal(normal) - get normals of each face in order;
//                      returns false when finished all faces;
//   GetNextUnitNormal(normal) - get normals of unit length of each face
//                      in order; returns false when finished all faces;
//   GetNextEdgeIndeces(i1, i2, edgeFlag) - get indeces of the next edge;
//                      returns false for the last edge;
//   GetNextEdge(p1, p2, edgeFlag) - get next edge;
//                      returns false for the last edge;
//   SetNumberOfRotationSteps(int n) - Set number of steps for whole circle;

// History:
// 21st February 2000  Evgeni Chernaev, John Allison
// - Re-written to inherit HepPolyhedron.
//
// 11.03.05 J.Allison
// - Added fNumberOfRotationStepsAtTimeOfCreation and access method.
//   (NumberOfRotationSteps is also called number of sides per circle or
//   line segments per circle - see
//   /vis/viewer/set/lineSegmentsPerCircle.)
// 20.06.05 G.Cosmo
// - Added GPolyhedronEllipsoid.
// 09.03.06 J.Allison
// - Added operator<<.

// #include "globals.hh"
#include "HepPolyhedron.h"
// #include "G4Visible.hh"

class GPolyhedron : public HepPolyhedron/*, public G4Visible*/ {
public:
  GPolyhedron ();
  GPolyhedron (const HepPolyhedron& from);
  // Use compiler defaults for copy contructor and assignment.  (They
  // invoke their counterparts in HepPolyhedron and G4Visible.)
  virtual ~GPolyhedron ();

  int GetNumberOfRotationStepsAtTimeOfCreation() const {
    return fNumberOfRotationStepsAtTimeOfCreation;
  }
private:
  int fNumberOfRotationStepsAtTimeOfCreation;
};

class GPolyhedronBox: public GPolyhedron {
public:
  GPolyhedronBox (double dx, double dy, double dz);
  virtual ~GPolyhedronBox ();
};

class GPolyhedronCone: public GPolyhedron {
public:
  GPolyhedronCone (double Rmn1, double Rmx1,
                    double Rmn2, double Rmx2, double Dz);
  virtual ~GPolyhedronCone ();
};

class GPolyhedronCons: public GPolyhedron {
public:
  GPolyhedronCons (double Rmn1, double Rmx1,
                    double Rmn2, double Rmx2, double Dz,
                    double Phi1, double Dphi);
  virtual ~GPolyhedronCons ();
};

class GPolyhedronPara: public GPolyhedron {
public:
  GPolyhedronPara (double Dx, double Dy, double Dz,
                    double Alpha, double Theta, double Phi);
  virtual ~GPolyhedronPara ();
};

class GPolyhedronPcon: public GPolyhedron {
public:
  GPolyhedronPcon (double phi, double dphi, int nz,
                    const double *z,
                    const double *rmin,
                    const double *rmax);
  GPolyhedronPcon (double phi, double dphi,
                    const std::vector<CLHEP::Hep2Vector> &rz);
  virtual ~GPolyhedronPcon ();
};

class GPolyhedronPgon: public GPolyhedron {
public:
  GPolyhedronPgon (double phi, double dphi, int npdv, int nz,
                    const double *z,
                    const double *rmin,
                    const double *rmax);
  GPolyhedronPgon (double phi, double dphi, int npdv,
                    const std::vector<CLHEP::Hep2Vector> &rz);

  virtual ~GPolyhedronPgon ();
};

class GPolyhedronSphere: public GPolyhedron {
public:
  GPolyhedronSphere (double rmin, double rmax,
                      double phi, double dphi,
                      double the, double dthe);
  virtual ~GPolyhedronSphere ();
};

class GPolyhedronTet: public GPolyhedron {
public:
  GPolyhedronTet (const double p0[3],
                   const double p1[3],
                   const double p2[3],
                   const double p3[3]);
  virtual ~GPolyhedronTet ();
};

class GPolyhedronTorus: public GPolyhedron {
public:
  GPolyhedronTorus (double rmin, double rmax, double rtor,
                    double phi, double dphi);
  virtual ~GPolyhedronTorus ();
};

class GPolyhedronTrap: public GPolyhedron {
public:
  GPolyhedronTrap (double Dz, double Theta, double Phi,
                    double Dy1,
                    double Dx1, double Dx2, double Alp1,
                    double Dy2,
                    double Dx3, double Dx4, double Alp2);
  virtual ~GPolyhedronTrap ();
};

class GPolyhedronTrd1: public GPolyhedron {
public:
  GPolyhedronTrd1 (double Dx1, double Dx2,
                    double Dy, double Dz);
  virtual ~GPolyhedronTrd1 ();
};

class GPolyhedronTrd2: public GPolyhedron {
public:
  GPolyhedronTrd2 (double Dx1, double Dx2,
                    double Dy1, double Dy2, double Dz);
  virtual ~GPolyhedronTrd2 ();
};

class GPolyhedronTube: public GPolyhedron {
public:
  GPolyhedronTube (double Rmin, double Rmax, double Dz);
  virtual ~GPolyhedronTube ();
};

class GPolyhedronTubs: public GPolyhedron {
public:
  GPolyhedronTubs (double Rmin, double Rmax, double Dz,
                    double Phi1, double Dphi);
  virtual ~GPolyhedronTubs ();
};

class GPolyhedronParaboloid: public GPolyhedron {
 public:
  GPolyhedronParaboloid(double r1, double r2, double dz,
                         double sPhi, double dPhi);
  virtual ~GPolyhedronParaboloid ();
};

class GPolyhedronHype: public GPolyhedron {
 public:
  GPolyhedronHype(double r1, double r2, double tan1,
                   double tan2, double halfZ);
  virtual ~GPolyhedronHype ();
};

class GPolyhedronEllipsoid : public GPolyhedron {
 public:
  GPolyhedronEllipsoid(double dx, double dy, double dz,
                        double zcut1, double zcut2);
  virtual ~GPolyhedronEllipsoid ();
};

class GPolyhedronEllipticalCone : public GPolyhedron {
 public:
  GPolyhedronEllipticalCone(double dx, double dy, double z,
                             double zcut1);
  virtual ~GPolyhedronEllipticalCone ();
};

class GPolyhedronHyperbolicMirror : public GPolyhedron {
 public:
  GPolyhedronHyperbolicMirror(double a, double h, double r);
  virtual ~GPolyhedronHyperbolicMirror ();
};

std::ostream& operator<<(std::ostream& os, const GPolyhedron&);

#endif /* GPolyhedron_HH */
