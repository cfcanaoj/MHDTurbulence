# Copyright Spack Project Developers. See COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)

from spack_repo.builtin.build_systems.makefile import MakefilePackage
from spack.package import *


class Mhdturbulence(MakefilePackage):
    """3D magneto-hydrodynamic decaying turbulence using MPI and Fortran."""

    homepage = "https://github.com/cfcanaoj/MHDTurbulence"
    git = "https://github.com/cfcanaoj/MHDTurbulence.git"

    license("GPL-3.0", checked_by="freifrauvonbleifrei")

    version("main", branch="main")

    depends_on("c", type="build")
    depends_on("cxx", type="build")
    depends_on("fortran", type=("build", "link"))
    depends_on("nvhpc+mpi", type=("build", "link"))

    # OpenACC source requires the NVIDIA HPC SDK compiler
    requires(
        "%nvhpc",
        msg="MHDTurbulence OpenACC code requires the NVIDIA HPC SDK compiler (nvhpc)",
    )
    def edit(self, spec, prefix):
        """Edit Makefile to use Spack's compiler and flags."""
        makefile = join_path(self.stage.source_path, "srcacc", "Makefile")
        f = FileFilter(makefile)

        # Replace the Fortran compiler by spack-supplied mpif90
        f.filter(r"^\s*fc\s*=.*", "fc = mpif90")

        # Replace compiler options to use Spack's flags
        if "%oneapi" not in spec:
            f.filter( # compiler flags
                r"^\s*fopt\s*=.*", "fopt = -g -O3 -mcmodel=medium"
            )  # TODO optimized flags -03 also for intel?
            f.filter( # OpenACC flags
                r"^\s*foptopenacc\s*=.*", "foptopenacc = -fopenacc"
            )
            # Replace Intel's OpenMP flag by most common one
            f.filter("-qopenmp", "-fopenmp")

    def build(self, spec, prefix):
        with working_dir(join_path(self.stage.source_path, "srcacc")):
            make()

    def install(self, spec, prefix):
        mkdirp(prefix.bin)
        install("exe/Simulation.x", prefix.bin)
