#!/bin/bash

rm *.vtu; rm *.vtk
./tp-imef.py TP-imef.msh
./tp-imef.py TP-imef-9nodos.msh