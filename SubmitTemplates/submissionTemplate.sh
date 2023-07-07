#!/bin/bash

export ACCGROUP="group_u_NA61.u_wj"

condor_submit /afs/cern.ch/user/c/camussol/private/Optics/SingleFile/%(subFileName)s.sub
