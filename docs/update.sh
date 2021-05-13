#!/bin/sh
rm -r ./_build
sphinx-apidoc -f -o ./source ../opentiva
make clean html

