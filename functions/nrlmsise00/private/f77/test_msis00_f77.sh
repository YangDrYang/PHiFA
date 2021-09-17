#!/bin/bash
clear

rm -f nrlmsise-test

ifort m_nrlmsise00.f nrlmsise00_driver.f -o nrlmsise-test

rm -f *.mod *.o

./nrlmsise-test