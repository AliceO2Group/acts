target_sources(
    ActsCore
    PRIVATE Common.cpp Direction.cpp ParticleData.cpp ParticleIdHelper.cpp
)

include(ActsCodegen)

if(ACTS_ADDITIONAL_PARTICLES_CSV)
message(STATUS "Using additional particles CSV: ${ACTS_ADDITIONAL_PARTICLES_CSV}")
acts_code_generation(
    ADD_TO_TARGET ActsCore
    PYTHON ${PROJECT_SOURCE_DIR}/Fatras/scripts/generate_particle_data_table.py
    PYTHON_ARGS --additional-particles ${ACTS_ADDITIONAL_PARTICLES_CSV}
    KEY Core/src/Definitions/ParticleDataTable.hpp
)
else()
acts_code_generation(
    ADD_TO_TARGET ActsCore
    KEY Core/src/Definitions/ParticleDataTable.hpp
)
endif()
