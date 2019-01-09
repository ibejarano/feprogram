#!/bin/bash

rm *.vtu; rm *.vtk
./main.py TP-imef.msh
./main.py TP-imef-9nodos.msh