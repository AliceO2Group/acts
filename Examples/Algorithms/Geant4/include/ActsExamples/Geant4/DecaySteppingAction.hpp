// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <memory>

#include <G4UserSteppingAction.hh>

class G4Step;

namespace ActsExamples::Geant4 {

class EventStore;

/// @class DecaySteppingAction
///
/// @brief Captures daughter particles at the exact moment of decay
///
/// This stepping action intercepts decay processes and stores information
/// about daughter particles at their creation point (decay vertex), before
/// any propagation through the detector material.
class DecaySteppingAction final : public G4UserSteppingAction {
 public:
  /// @brief Configuration of the DecaySteppingAction
  struct Config {
    /// Event store to write to
    EventStore* eventStore = nullptr;
  };

  /// Construct the stepping action
  ///
  /// @param cfg the configuration struct
  /// @param logger the ACTS logging instance
  DecaySteppingAction(const Config& cfg,
                      std::unique_ptr<const Acts::Logger> logger =
                          Acts::getDefaultLogger("DecaySteppingAction",
                                                  Acts::Logging::INFO));
  ~DecaySteppingAction() override = default;

  /// Action per step
  ///
  /// @param step the Geant4 step of the particle
  void UserSteppingAction(const G4Step* step) override;

 protected:
  /// Access the event store
  EventStore& eventStore() const { return *m_cfg.eventStore; }

  /// The logger instance
  const Acts::Logger& logger() const { return *m_logger; }

 private:
  /// The configuration
  Config m_cfg;
  /// The logging instance
  std::unique_ptr<const Acts::Logger> m_logger;
};

}  // namespace ActsExamples::Geant4
