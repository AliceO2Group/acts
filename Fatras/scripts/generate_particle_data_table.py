#!/usr/bin/env python3

# /// script
# dependencies = [
#   "particle==0.24.0",
# ]
# ///

#
# use scikit-hep/particle to generate c++ code for the particle data table.
#

import io
import sys
import subprocess
import argparse
from pathlib import Path

from particle import Particle


CODE_HEADER = """\
// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// The entries within this file have been automatically created using the
// particle data files from Review of Particle Physics
// by the Berkeley Particle Data Group.

#pragma once

#include <cstdint>
#include <map>
#include <limits>

"""

# Default list of additional particles, e.g. BSM particles that are not in the
# PDG tables of scikit-hep/particle. Their PDG IDs and names are compiled in;
# their masses and charges can be overridden at runtime, see CSV_FILE_READER.
DEFAULT_ADDITIONAL_PARTICLES_CSV = Path(__file__).parent / "additional_particles.csv"

# Only emitted if there are additional particles. The masses and charges of the
# additional particles are compiled in as defaults, and can be replaced at
# runtime (e.g. to scan a mass range without rebuilding) by a CSV file found as
#   1. $ACTS_ADDITIONAL_PARTICLES_CSV, if set (it is an error if unreadable)
#   2. $ACTS_ROOT/share/Acts/additional_particles.csv, if it exists
# The file has the same format as the build-time CSV. Particles are matched by
# PDG ID, so it only needs to contain the particles whose values change.
CSV_FILE_READER = """
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct AdditionalParticleValues {
  float charge;  // in e
  float mass;    // in GeV
};

inline std::map<std::int32_t, AdditionalParticleValues> readAdditionalParticlesCsv(
    const std::string& csvPath,
    const std::map<std::int32_t, AdditionalParticleValues>& known) {
  std::ifstream input(csvPath);
  if (!input.is_open()) {
    throw std::runtime_error("Failed to open additional particles CSV file " + csvPath);
  }

  std::map<std::int32_t, AdditionalParticleValues> values;
  std::string line;
  for (int lineNumber = 1; std::getline(input, line); ++lineNumber) {
    const auto first = line.find_first_not_of(" \\t\\r");
    if (first == std::string::npos || line[first] == '#') {
      continue;  // skip empty lines and comments
    }

    std::vector<std::string> columns;
    std::stringstream lineStream(line);
    std::string value;
    while (std::getline(lineStream, value, ',')) {
      columns.push_back(value);
    }

    const std::string where = " at line " + std::to_string(lineNumber) + " in file " + csvPath;
    if (columns.size() != 4) {
      throw std::runtime_error("Expected 4 columns (ID,Name,Mass,Charge)" + where);
    }

    std::int32_t pdg = 0;
    AdditionalParticleValues particle{};
    try {
      pdg = std::stoi(columns[0]);
      particle.mass = std::stof(columns[2]) / 1000.f;    // MeV -> GeV
      particle.charge = std::stof(columns[3]) / 3.f;     // three-charge -> e
    } catch (...) {
      throw std::runtime_error("Failed to convert ID, mass or charge to a number" + where);
    }
    if (!known.contains(pdg)) {
      throw std::runtime_error(
          "PDG ID " + columns[0] + where +
          " is not one of the additional particles ACTS was built with; "
          "add it to Fatras/scripts/additional_particles.csv and rebuild");
    }
    values[pdg] = particle;
  }
  return values;
}

"""

# Emitted after the compiled-in defaults (additionalParticlesDefaults).
ADDITIONAL_PARTICLES_LOOKUP = """
inline const std::map<std::int32_t, AdditionalParticleValues>& additionalParticles() {
  static const std::map<std::int32_t, AdditionalParticleValues> values = [] {
    std::map<std::int32_t, AdditionalParticleValues> result = additionalParticlesDefaults();

    std::string csvPath;
    if (const char* env = std::getenv("ACTS_ADDITIONAL_PARTICLES_CSV"); env != nullptr && *env != '\\0') {
      csvPath = env;
    } else if (const char* root = std::getenv("ACTS_ROOT"); root != nullptr && *root != '\\0') {
      const std::string installed = std::string(root) + "/share/Acts/additional_particles.csv";
      if (std::ifstream(installed).good()) {
        csvPath = installed;
      }
    }

    if (!csvPath.empty()) {
      for (const auto& [pdg, particle] : readAdditionalParticlesCsv(csvPath, result)) {
        result[pdg] = particle;
      }
    }
    return result;
  }();
  return values;
}

"""


def read_additional_particles(csv_path):
    """
    Read the additional particles CSV. Returns a list of
    ``(pdgid, name, mass_MeV, three_charge)`` tuples.
    """
    particles = []
    with open(csv_path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue  # skip empty lines and comments
            parts = line.split(",")
            if len(parts) != 4:
                print("Warning: Skipping invalid line in additional particles CSV:", line, file=sys.stderr)
                continue
            pdgid_str, name, mass_str, charge_str = parts
            try:
                particles.append((int(pdgid_str), name, float(mass_str), float(charge_str)))
            except ValueError:
                print(
                    "Warning: Skipping line with invalid number format in additional particles CSV: ", line,
                    file=sys.stderr)
    return particles


def convert_mass_to_GeV(mass_MeV):
    """
    Convert the mass from MeV to GeV.
    """
    # return mass_MeV
    return mass_MeV / 1000.0


def convert_charge_to_e(charge):
    """
    Convert the charge from three-charge to electron charge.
    """
    return charge / 3.0


def identity(v):
    """
    Return the value unchanged.
    """
    return v


def generate_code(additional_particles_csv: Path | None = DEFAULT_ADDITIONAL_PARTICLES_CSV):
    """
    Generate C++ source text containing particle lookup maps.

    By default, the generated code contains all particles provided by
    ``particle.Particle.all()`` and creates one ``std::map`` per property:
    charge, mass, and name, keyed by signed PDG ID.

    If ``additional_particles_csv`` exists, the function also appends the
    particles from that CSV to each generated map. Their names are fixed at
    build time; their masses and charges are compiled in as defaults that can
    be overridden at runtime, see ``CSV_FILE_READER``.

    Parameters
    ----------
    additional_particles_csv : str | None, optional
        Path to a CSV file containing additional particles. Defaults to
        ``additional_particles.csv`` next to this script.

        Expected row format (no header required):
            ``ID,Name,Mass,Charge``

        Where:
        - ``ID`` is a signed integer PDG ID.
        - ``Name`` is a particle name string.
        - ``Mass`` is in MeV.
        - ``Charge`` is in three-charge units.

        Empty lines and lines starting with ``#`` are ignored.
        Invalid rows are skipped with a warning to stderr.

    Returns
    -------
    str
        Generated C++ code as a single string, ending with a newline.
    """

    lines = CODE_HEADER.splitlines()

    additional_particles = []
    if additional_particles_csv is not None and Path(additional_particles_csv).exists():
        additional_particles = read_additional_particles(additional_particles_csv)

    if additional_particles:
        lines.append(CSV_FILE_READER)
        lines.append(
            "inline const std::map<std::int32_t, AdditionalParticleValues>& additionalParticlesDefaults() {"
        )
        lines.append(
            "static const std::map<std::int32_t, AdditionalParticleValues> map = {"
        )
        for pdgid, name, mass, three_charge in additional_particles:
            lines.append(f"// {name}")
            lines.append(
                f"{{ {pdgid}, {{ {convert_charge_to_e(three_charge)}f, {convert_mass_to_GeV(mass)}f }} }},"
            )
        lines.append("};")
        lines.append("return map;")
        lines.append("}")
        lines.append(ADDITIONAL_PARTICLES_LOOKUP)

    lines.append(
        f"static constexpr std::uint32_t kParticlesCount = {len(Particle.all())}u;"
    )

    # ensure the rows are sorted by the signed pdg number (first column)
    # name, c++ type, and output format for each column
    columns = [
        ("three_charge", "Charge", "float", "{}", convert_charge_to_e),
        ("mass", "Mass", "float", "{}f", convert_mass_to_GeV),
        ("name", "Name", " const  char* const", '"{}"', identity),
    ]
    # build a separate array for each column
    for variable_name, target_name, type_name, value_format, transform in columns:
        cpp_name = f"particlesMap{target_name}"
        map_type = f"std::map<std::int32_t, {type_name}>"

        lines.append(f"inline const {map_type}& {cpp_name}() {{")
        lines.append(f"static const {map_type} map = {{")

        for p in Particle.all():
            value = getattr(p, variable_name)

            lines.append(f"// {p.name}")

            if value is None and type_name == "float":
                lines.append(
                    f"{{ {int(p.pdgid)} , std::numeric_limits<float>::quiet_NaN() }},"
                )
            else:
                lines.append(
                    f"{{ {int(p.pdgid)}, {value_format.format(transform(value))} }},"
                )

        # Masses and charges of additional particles are looked up at runtime,
        # names are fixed at build time
        for pdgid, name, _, _ in additional_particles:
            lines.append(f"// {name}")
            if variable_name == "name":
                lines.append(f"{{ {pdgid}, {value_format.format(name)} }},")
            elif variable_name == "mass":
                lines.append(f"{{ {pdgid}, additionalParticles().at({pdgid}).mass }},")
            elif variable_name == "three_charge":
                lines.append(f"{{ {pdgid}, additionalParticles().at({pdgid}).charge }},")
            else:
                raise ValueError(f"Unsupported variable name: {variable_name}")

        lines.append("};")
        lines.append("return map;")
        lines.append("}")

    # ensure we end with a newline
    lines.append("")
    return "\n".join(lines)


def clang_format(content):
    """
    Format the given content using clang-format and return it.
    """
    args = [
        "clang-format",
        "--assume-filename=ParticleData.hpp",
        "-",
    ]
    process = subprocess.run(
        args,
        input=content,
        capture_output=True,
        check=True,
        encoding="utf-8",
        text=True,
    )
    return process.stdout


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Generate the particle data table.")
    p.add_argument("output", type=Path, nargs="?", default=None, help="Output file.")
    p.add_argument(
        "--format", action="store_true", help="Run clang-format on the output."
    )
    p.add_argument(
        "--additional-particles",
        type=Path,
        help="CSV file with additional particles.",
        default=DEFAULT_ADDITIONAL_PARTICLES_CSV,
    )

    args = p.parse_args()
    if args.output is None:
        output_file = sys.stdout
    else:
        # will overwrite existing file
        output_file = io.open(args.output, mode="wt", encoding="utf-8")

    code = generate_code(additional_particles_csv=args.additional_particles)
    if args.format:
        code = clang_format(code)
    output_file.write(code)
