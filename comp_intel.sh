#!/bin/bash

#module purge

module load intel intel-mpi boost

make $1
