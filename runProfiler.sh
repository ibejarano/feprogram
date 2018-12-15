#!/bin/bash

rm *.vtu; rm *.vtk
kernprof -l -v ./tp-imef.py TP-imef-9nodos.msh
