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
/// \file B4/B4a/src/PrimaryGeneratorAction.cc
/// \brief Implementation of the B4::PrimaryGeneratorAction class
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4AnalysisManager.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include <fstream>
#include <sstream>
#include <vector>
#include <stdexcept>

namespace B4
{
    // Vector global (compartido por todas las instancias) para guardar energías
    static std::vector<G4double> energias;

    PrimaryGeneratorAction::PrimaryGeneratorAction()
    {
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);

        auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        fParticleGun->SetParticleDefinition(particleDefinition);
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

        // Cargar energías solo una vez
        /*
        if (energias.empty()) {
            std::ifstream infile("1e6.txt");
            if (!infile.is_open()) {
                G4Exception("PrimaryGeneratorAction", "FileError", FatalException, "No se pudo abrir 1e6.txt");
            }

            std::string line;
            while (std::getline(infile, line)) {
                if (!line.empty()) {
                    energias.push_back(std::stod(line));
                }
            }

            if (energias.empty()) {
                G4Exception("PrimaryGeneratorAction", "EmptyFile", FatalException, "Archivo 1e6.txt vacío.");
            }
        }
        */
    }

    PrimaryGeneratorAction::~PrimaryGeneratorAction()
    {
        delete fParticleGun;
    }

    void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
    {
        G4int eventID = anEvent->GetEventID();
        G4double energyMeV = energias[eventID % energias.size()];

        auto analysisManager = G4AnalysisManager::Instance();
        //analysisManager->FillH1(2, energyMeV); // Energía en MeV
        //analysisManager->FillNtupleDColumn(2, energyMeV);
        //analysisManager->AddNtupleRow();

		// Set gun energy
        //fParticleGun->SetParticleEnergy(energyMeV * MeV);
        fParticleGun->SetParticleEnergy(1 * MeV);

        // Set gun position
        G4double xmin = -5 * cm;
        G4double xmax = +5 * cm;
        G4double randomx = xmin + (xmax - xmin) * G4UniformRand();

        G4double ymin = -5 * cm;
        G4double ymax = +5 * cm;
        G4double randomy = ymin + (ymax - ymin) * G4UniformRand();

		G4double zpos = 0.5 * cm;
        G4double xpos = 10 * cm;
        //fParticleGun->SetParticlePosition(G4ThreeVector(randomx, randomy, 0)); // Random inside detector
		//fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, 0));             // Centered inside detector
		fParticleGun->SetParticlePosition(G4ThreeVector(0, 0, zpos));       // Pointing to detector

        // Set gun direction
        G4double detsize = 5 * cm;
		G4double ratio = xpos / std::sqrt(detsize*detsize + xpos*xpos);
        //G4double cosTheta = (1 - ratio) * G4UniformRand() + ratio, phi = twopi * G4UniformRand(); 
        G4double cosTheta = 2 * G4UniformRand() - 1, phi = twopi * G4UniformRand();
        G4double sinTheta = std::sqrt(1. - cosTheta * cosTheta);
        G4double uz = sinTheta * std::cos(phi),
            uy = sinTheta * std::sin(phi),
            ux = cosTheta;
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, -1.0)); // Pointing to detector
        //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-ux, uy, uz));   // Solid angle
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

}

/*
namespace B4
{
    PrimaryGeneratorAction::PrimaryGeneratorAction()
    {
        G4int nofParticles = 1;
        fParticleGun = new G4ParticleGun(nofParticles);

        auto particleDefinition = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
        fParticleGun->SetParticleDefinition(particleDefinition);

        // Dirección fija hacia +Z (puedes cambiarla abajo si quieres aleatoriedad)
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    }

    PrimaryGeneratorAction::~PrimaryGeneratorAction()
    {
        delete fParticleGun;
    }

    void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
    {
        // Energía aleatoria entre 1 eV y 10 eV
        G4double emin = 10000 * eV;
        G4double emax = 100000 * eV;
        G4double randomE = emin + (emax - emin) * G4UniformRand();

        fParticleGun->SetParticleEnergy(randomE);

        // Posición aleatoria en un cuadrado de 30x60 cm
        G4double randomx = (-15 * cm) + 30 * cm * G4UniformRand();
        G4double randomy = (-30 * cm) + 60 * cm * G4UniformRand();
        G4double posZ = -12. * m;

        fParticleGun->SetParticlePosition(G4ThreeVector(randomx, randomy, posZ));

        // Dirección hacia +Z (puedes hacerla hemisférica aleatoria si quieres)
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

        fParticleGun->GeneratePrimaryVertex(anEvent);
    }

}
*/