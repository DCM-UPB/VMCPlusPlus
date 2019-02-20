[![Build Status](https://travis-ci.com/DCM-UPB/VMCPlusPlus.svg?branch=master)](https://travis-ci.com/DCM-UPB/VMCPlusPlus)
[![Coverage Status](https://coveralls.io/repos/github/DCM-UPB/VMCPlusPlus/badge.svg?branch=master)](https://coveralls.io/github/DCM-UPB/VMCPlusPlus?branch=master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/e4e5ff4974e24a4faffff3d6c7010852)](https://www.codacy.com/app/NNVMC/VMCPlusPlus?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=DCM-UPB/VMCPlusPlus&amp;utm_campaign=Badge_Grade)
[![CodeFactor](https://www.codefactor.io/repository/github/dcm-upb/vmcplusplus/badge)](https://www.codefactor.io/repository/github/dcm-upb/vmcplusplus)

# VariationalMonteCarloPlusPlus

C++ Library for doing a Variational Monte Carlo simulation.
It is built upon the MCI++ (https://github.com/DCM-UPB/MCIntegratorPlusPlus) and NoisyFunMin (https://github.com/DCM-UPB/NoisyFunMin) libraries.
It also uses the GSL (GNU Scientific Library).

In `doc/` there is a user manual in pdf.

In `examples/` there are several examples.



# Build the library

We use the CMake build system, so you need to have it on your system to build the library out of the box.
Then copy the file `config_template.sh` to `config.sh`, edit it to your liking and then simply execute the command

   `./build.sh`

Note that we build out-of-tree, so the compiled library and executable files can be found in the directories under `./build/`.
