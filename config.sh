#!/bin/bash
mkdir ./input
scp wasa:/data7/users/khreptak/OUTPUT/DATA/DATA-newcuts-AddGammaCut-offset-bound-pdpi0.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/MC-newcuts-AddGammaCut-pd-bound-pdpi0.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/MC-newcuts-AddGammaCut-pd-pdpi0.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/MC-newcuts-AddGammaCut-x6-pd-bound-pdpi0.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/MC-newcuts-AddGammaCut-x6-pd-pdpi0.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/MC-acceptance.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/MC-newcuts-AddGammaCut-pd-bound-pdpi0_G10.root ./input
scp wasa:/data7/users/khreptak/OUTPUT/MC/LUMIN/MC-dnpi_plus_thetacut-x6.root ./input
mkdir ./output
mkdir ./output/plots
