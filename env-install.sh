export G4BUILD_DIR=/home/zack/Desktop/geant4.9.6.p02-installopt/share/Geant4-9.6.2/geant4make
export G4INSTALL_DIR=/home/zack/Desktop/geant4.9.6.p02-installopt
export G4WORKDIR=$G4INSTALL_DIR/Gmake-Workdir	
export G4SYSTEM=Linux-g++
export G4INSTALL=/home/zack/Desktop/geant4.9.6.p02-installopt

export G4BIN=$G4WORKDIR/bin/Shooter

if [ -x $G4BUILD_DIR/geant4make.sh ]; then
   source $G4BUILD_DIR/geant4make.sh
else
   echo "Cannot find shell config file geant4make.csh in directory $G4BUILD_DIR "
   exit 1
fi
if [ -x  $G4INSTALL_DIR/bin/geant4.sh ]; then
    source $G4INSTALL_DIR/bin/geant4.sh
else
   echo "Cannot find shell config file geant4.csh in directory $G4INSTALL_DIR "
   exit 2
fi
