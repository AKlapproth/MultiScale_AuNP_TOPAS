## General information

0) Pre-requisites:

TOPAS installed with recommended `OS system, c++ and cmake versions`, see 
   topas https://topas.readthedocs.io/en/latest/getting-started/install.htm

1) We recommend having a global directory for extensions named topas_extension and move in TOPAS-nBio there 

   Linux:
	mkdir ~/topas_extensions

   Mac: 
        mkdir /Applications/topas_extensions

2) Unzip TOPAS-nBio directory in topas_extensions and navigate to the topas directory

  Linux:
        cd ~/topas
  Mac:
	cd /Applications/topas

3) Unzip the Geant4Headers.zip

  Linux:
	unzip -e Geant4Headers.zip
  Mac:
	unzip -e Geant4Headers.zip
        
4) Build the extensions

  Linux:
	cmake ./ -DTOPAS_EXTENSIONS_DIR=~/topas_extensions/TOPAS-nBio
	make -j4
  Mac:
	cmake ./ -DTOPAS_EXTENSIONS_DIR=/Applications/topas_extensions/TOPAS-nBio
	make -j4
 
5) Run the demos. For some demos, a pause before quit is enabled, then, write `exit` at the terminal prompt.

  Linux:
	source rundemos.csh

  Mac:
	source rundemos.csh

