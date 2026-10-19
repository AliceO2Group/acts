target_sources(
    ActsCore
    PRIVATE Common.cpp Direction.cpp ParticleData.cpp ParticleIdHelper.cpp
)

include(ActsCodegen)

acts_code_generation(
    ADD_TO_TARGET ActsCore
    KEY Core/src/Definitions/ParticleDataTable.hpp
)

# Installed as the runtime default for the masses and charges of the additional
# particles, see Fatras/scripts/generate_particle_data_table.py
install(
    FILES ${PROJECT_SOURCE_DIR}/Fatras/scripts/additional_particles.csv
    DESTINATION ${CMAKE_INSTALL_DATADIR}/Acts
)
