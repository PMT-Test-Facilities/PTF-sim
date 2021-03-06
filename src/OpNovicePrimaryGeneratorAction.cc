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
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "OpNovicePrimaryGeneratorAction.hh"
#include "OpNovicePrimaryGeneratorMessenger.hh"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNovicePrimaryGeneratorAction::OpNovicePrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(), 
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
  fParticleGPS = new G4GeneralParticleSource();

  fSourceType = "gun";
  
  //create a messenger for this class
  fGunMessenger = new OpNovicePrimaryGeneratorMessenger(this);

  //default kinematic
  //
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("e+");

  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleTime(0.0*ns);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0*cm,0.0*cm,0.0*cm));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
  fParticleGun->SetParticleEnergy(500.0*keV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

OpNovicePrimaryGeneratorAction::~OpNovicePrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fGunMessenger;
  delete fParticleGPS;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNovicePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  static int scan_location = 0;

  const double scan_xmin = -250.0 * CLHEP::mm;
  const double scan_xmax =  250.0 * CLHEP::mm;
  const double scan_ymin = -250.0 * CLHEP::mm;
  const double scan_ymax =  250.0 * CLHEP::mm;
  const int scan_nx = 250;
  const int scan_ny = 250;
  const int scan_N = scan_nx * scan_ny;

  const double scan_dx = (scan_xmax - scan_xmin) / scan_nx;
  const double scan_dy = (scan_ymax - scan_ymin) / scan_ny;

  if ( fSourceType == "gun" ){
    fParticleGun->GeneratePrimaryVertex(anEvent);

  } else if ( fSourceType == "vert-inc-scan" ){
    
    fParticleGun->SetNumberOfParticles( 10000 );
    fParticleGun->SetParticleDefinition( G4ParticleTable::GetParticleTable()->FindParticle("opticalphoton") );
    int nparticles = fParticleGun->GetNumberOfParticles();
    //std::cout<<"particles="<<fParticleGun->GetParticleDefinition()->GetParticleName()<<std::endl;
    G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
    //std::cout<<"direction="<<dir.x()<<" , "<<dir.y()<<", "<<dir.z()<<std::endl;

    int ix = scan_location / scan_nx;
    int iy = scan_location % scan_nx;

    double xloc = scan_xmin + ix * scan_dx;
    double yloc = scan_ymin + iy * scan_dy;
    
    G4ThreeVector loc = fParticleGun->GetParticlePosition();

    fParticleGun->SetParticlePosition( G4ThreeVector( xloc, yloc, loc.z() ) );
    fParticleGun->GeneratePrimaryVertex(anEvent);
    scan_location = (++scan_location) % scan_N;
    //if ( scan_location % 100 == 0 ){
      std::cout<<"scan_location="<<scan_location<<" of "<<scan_N<<std::endl;
      std::cout<<"Generating "<<nparticles<<" at "<<xloc/CLHEP::cm<<", "<<yloc/CLHEP::cm<<", "<<loc.z()/CLHEP::cm<<"  cm"<<std::endl;
      //}

  } else if ( fSourceType == "norm-inc-scan" ){


  } else if ( fSourceType == "angle-scan" ){
      
  } else if ( fSourceType == "gps" ) {
    std::cout<<"OpNovicePrimaryGeneratorAction: using gps source"<<std::endl;
    int numpart = fParticleGPS->GetNumberOfParticles();
    std::cout<<"Generate "<<numpart<<" particles"<<std::endl;
    fParticleGPS->GeneratePrimaryVertex(anEvent);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNovicePrimaryGeneratorAction::SetOptPhotonPolar()
{
 G4double angle = G4UniformRand() * 360.0*deg;
 SetOptPhotonPolar(angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void OpNovicePrimaryGeneratorAction::SetOptPhotonPolar(G4double angle)
{
 if (fParticleGun->GetParticleDefinition()->GetParticleName()!="opticalphoton")
   {
     G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
               "the particleGun is not an opticalphoton" << G4endl;
     return;
   }

 G4ThreeVector normal (1., 0., 0.);
 G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
 G4ThreeVector product = normal.cross(kphoton);
 G4double modul2       = product*product;
 
 G4ThreeVector e_perpend (0., 0., 1.);
 if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
 G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
 
 G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
 fParticleGun->SetParticlePolarization(polar);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......