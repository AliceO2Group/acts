// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Trajectories.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsFatras/EventData/Particle.hpp"

#include <functional>

#include "edm4hep/MCParticle.h"
#include "edm4hep/MutableMCParticle.h"
#include "edm4hep/MutableSimTrackerHit.h"
#include "edm4hep/MutableTrack.h"
#include "edm4hep/MutableTrackerHit.h"
#include "edm4hep/MutableTrackerHitPlane.h"
#include "edm4hep/SimTrackerHit.h"
#include "edm4hep/TrackerHit.h"
#include "edm4hep/TrackerHitCollection.h"
#include "edm4hep/TrackerHitPlane.h"

namespace ActsExamples {
namespace EDM4hepUtil {

static constexpr std::int32_t EDM4HEP_ACTS_POSITION_TYPE = 42;

using MapParticleIdFrom =
    std::function<ActsFatras::Barcode(edm4hep::MCParticle particle)>;
using MapParticleIdTo =
    std::function<edm4hep::MCParticle(ActsFatras::Barcode particleId)>;

using MapGeometryIdFrom =
    std::function<Acts::GeometryIdentifier(std::uint64_t cellId)>;
using MapGeometryIdTo =
    std::function<std::uint64_t(Acts::GeometryIdentifier geometryId)>;

/// Reads a Fatras particle from EDM4hep.
///
/// Inpersistent information:
/// - particle ID
/// - process
ActsFatras::Particle readParticle(const edm4hep::MCParticle& from,
                                  const MapParticleIdFrom& particleMapper);

/// Write a Fatras particle into EDM4hep.
///
/// Inpersistent information:
/// - particle ID
/// - process
void writeParticle(const ActsFatras::Particle& from,
                   edm4hep::MutableMCParticle to);

/// Reads a Fatras hit from EDM4hep.
///
/// Inpersistent information:
/// - after4 momentum
/// - hit index
/// - digitization channel
ActsFatras::Hit readSimHit(const edm4hep::SimTrackerHit& from,
                           const MapParticleIdFrom& particleMapper,
                           const MapGeometryIdFrom& geometryMapper);

/// Writes a Fatras hit to EDM4hep.
///
/// Inpersistent information:
/// - after4 momentum
/// - hit index
/// - digitization channel
void writeSimHit(const ActsFatras::Hit& from, edm4hep::MutableSimTrackerHit to,
                 const MapParticleIdTo& particleMapper,
                 const MapGeometryIdTo& geometryMapper);

/// Reads a measurement cluster from EDM4hep.
///
/// Inpersistent information:
/// - hit index
/// - 1D local coords?
/// - segment path
///
/// Known issues:
/// - cluster channels are read from inappropriate fields
/// - local 2D coordinates and time are read from position
Measurement readMeasurement(const edm4hep::TrackerHitPlane& from,
                            const edm4hep::TrackerHitCollection* fromClusters,
                            Cluster* toCluster,
                            const MapGeometryIdFrom& geometryMapper);

/// Writes a measurement cluster to EDM4hep.
///
/// Inpersistent information:
/// - hit index
/// - 1D local coords?
/// - segment path
///
/// Known issues:
/// - cluster channels are written to inappropriate fields
/// - local 2D coordinates and time are written to position
void writeMeasurement(const Measurement& from,
                      edm4hep::MutableTrackerHitPlane to,
                      const Cluster* fromCluster,
                      edm4hep::TrackerHitCollection& toClusters,
                      const MapGeometryIdTo& geometryMapper);

/// Writes a trajectory to EDM4hep.
///
/// Inpersistent information:
/// - trajectory state incomplete
/// - relation to the particles
///
/// Known issues:
/// - curvature parameter
/// - track state local coordinates are written to (D0,Z0)
/// - covariance incorrect
void writeTrajectory(const Trajectories& from, edm4hep::MutableTrack to,
                     std::size_t fromIndex,
                     const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap);

}  // namespace EDM4hepUtil
}  // namespace ActsExamples
