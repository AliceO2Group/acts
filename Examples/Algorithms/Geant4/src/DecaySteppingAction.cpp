// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/DecaySteppingAction.hpp"

#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <G4ParticleDefinition.hh>
#include <G4Step.hh>
#include <G4Track.hh>
#include <G4VProcess.hh>

namespace ActsExamples::Geant4 {

DecaySteppingAction::DecaySteppingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserSteppingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void DecaySteppingAction::UserSteppingAction(const G4Step* step) {
  // Get the track and process information
  G4Track* track = step->GetTrack();
  const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
  
  if (process == nullptr) {
    return;
  }

  // Check if this is a decay process
  G4String processName = process->GetProcessName();
  if (processName != "Decay") {
    return;
  }

  // Get the parent particle barcode
  G4int trackId = track->GetTrackID();
  if (!eventStore().trackIdMapping.contains(trackId)) {
    return;
  }
  
  SimBarcode parentBarcode = eventStore().trackIdMapping.at(trackId);

  // Get the number of secondaries created in this step
  const G4TrackVector* secondaries = step->GetSecondary();
  if (secondaries == nullptr || secondaries->size() == 0) {
    return;
  }

  // Unit conversions G4->ACTS
  constexpr double convertTime = Acts::UnitConstants::ns / CLHEP::ns;
  constexpr double convertLength = Acts::UnitConstants::mm / CLHEP::mm;
  constexpr double convertEnergy = Acts::UnitConstants::GeV / CLHEP::GeV;

  // Get decay vertex information from the post-step point (where decay happens)
  G4ThreeVector decayPosition = convertLength * step->GetPostStepPoint()->GetPosition();
  G4double decayTime = convertTime * step->GetPostStepPoint()->GetGlobalTime();

  ACTS_VERBOSE("Decay detected: parent barcode " << parentBarcode 
               << ", parent trackID " << trackId
               << ", PDG " << track->GetParticleDefinition()->GetPDGEncoding()
               << ", " << secondaries->size() << " daughters"
               << " at position (" << decayPosition.x() << ", " 
               << decayPosition.y() << ", " << decayPosition.z() << ")");

  // Create a list to store decay daughters for this parent
  std::vector<DecayVertexInfo> decayDaughters;

  // Store information about each daughter particle at the decay vertex
  int daughterIndex = 0;
  for (const auto* secondary : *secondaries) {
    // Get daughter particle information
    const G4ParticleDefinition* particleDef = secondary->GetParticleDefinition();
    G4int pdg = particleDef->GetPDGEncoding();
    G4double charge = particleDef->GetPDGCharge();
    G4double mass = convertEnergy * particleDef->GetPDGMass();
    
    // Get momentum at creation (decay vertex)
    G4ThreeVector momentum = secondary->GetMomentum();
    G4double px = convertEnergy * momentum.x();
    G4double py = convertEnergy * momentum.y();
    G4double pz = convertEnergy * momentum.z();
    G4double p = convertEnergy * momentum.mag();

    if(p<1.0e-04)
    {
      ACTS_WARNING("Daughter " << daughterIndex << " has very low momentum p = " << p << " GeV. Skipping.");
      continue;
    }
    
    // Get direction
    G4ThreeVector direction = momentum.unit();

    // Create unique barcode for each daughter (will be used by ParticleTrackingAction)
    SimBarcode daughterBarcode = parentBarcode.makeDescendant();
    
    // Generate unique subparticle ID for this daughter
    auto key = EventStore::BarcodeWithoutSubparticle::Zeros();
    key.set(0, daughterBarcode.vertexPrimary())
        .set(1, daughterBarcode.vertexSecondary())
        .set(2, daughterBarcode.particle())
        .set(3, daughterBarcode.generation());
    daughterBarcode.setSubParticle(++eventStore().subparticleMap[key]);
    
    // Store decay vertex info for this daughter
    DecayVertexInfo decayInfo;
    decayInfo.barcode = daughterBarcode;
    decayInfo.parentBarcode = parentBarcode;  // Store parent barcode
    decayInfo.position = Acts::Vector3(decayPosition.x(), decayPosition.y(), decayPosition.z());
    decayInfo.time = decayTime;
    decayInfo.momentum = Acts::Vector3(px, py, pz);
    decayInfo.pdg = pdg;
    decayDaughters.push_back(decayInfo);
    
    // Create SimParticleState for the daughter at decay vertex for particlesDecay collection
    SimParticleState daughterState(daughterBarcode, Acts::PdgParticle{pdg}, charge, mass);
    daughterState.setPosition4(decayPosition.x(), decayPosition.y(), decayPosition.z(), decayTime);
    daughterState.setDirection(direction.x(), direction.y(), direction.z());
    daughterState.setAbsoluteMomentum(p);
    
    // Create a SimParticle with both initial and final state set to decay vertex
    SimParticle daughterParticle(daughterState, daughterState);
    
    // Store in particlesDecay collection
    eventStore().particlesDecay.insert(daughterParticle);
    
    // Store daughter-to-mother barcode mapping
    eventStore().daughterToMotherMap[daughterBarcode] = parentBarcode;
    
    ACTS_VERBOSE("  Daughter " << daughterIndex << ": barcode " << daughterBarcode 
                 << ", parent barcode " << parentBarcode
                 << ", PDG " << pdg 
                 << ", p = " << p << " GeV"
                 << ", px = " << px << ", py = " << py << ", pz = " << pz);
    
    daughterIndex++;
  }
  
  // Store the decay daughters indexed by parent track ID
  // ParticleTrackingAction will retrieve them when daughters start tracking
  eventStore().decayVertexMap[trackId] = decayDaughters;
}

}  // namespace ActsExamples::Geant4
