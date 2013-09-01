# $Id: GNUmakefile 69699 2013-05-13 08:50:30Z gcosmo $
# ----------------------------------------------------------------
# Makes test program in environment variable G4TARGET.
# ----------------------------------------------------------------
TESTTARGET=shooter;
ifndef G4TARGET
  G4TARGET := $(TESTTARGET)
endif

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

G4EXEC_BUILD = true

include $(G4INSTALL)/config/architecture.gmk

# Override some variables from binmake.gmk.
#
INCFLAGS := -I$(G4BASE)/geometry/management/include \
            -I$(G4BASE)/geometry/volumes/include \
            -I$(G4BASE)/geometry/navigation/include \
            -I$(G4BASE)/geometry/magneticfield/include \
            -I$(G4BASE)/geometry/solids/CSG/include \
            -I$(G4BASE)/global/management/include \
            -I$(G4BASE)/global/HEPRandom/include \
            -I$(G4BASE)/global/HEPGeometry/include \
	    -I$(G4BASE)/materials/include \
            -I$(G4BASE)/graphics_reps/include \
	    -I$(G4BASE)/intercoms/include

GLOBAL_LDLIBS   := -lG4geometry \
            -lG4graphics_reps \
            -lG4intercoms \
            -lG4global  

GRANULAR_LDLIBS   := \
	    -lG4csg \
            -lG4navigation \
	    -lG4volumes \
	    -lG4magneticfield \
	    -lG4geometrymng \
	    -lG4materials \
            -lG4graphics_reps \
	    -lG4intercoms \
            -lG4globman

ifdef G4GEOMETRY_VERBOSE
  CPPFLAGS += -DG4GEOMETRY_VERBOSE
endif

LDLIBS := $(GLOBAL_LDLIBS)
# LDLIBS := $(GRANULAR_LDLIBS)

ifdef G4LIB_USE_CLHEP
    LDLIBS += -lG4clhep
endif

#echobase:
#	echo "G4BASE is " $(G4BASE)
#echoinstall: 
#	echo "G4INSTALL is " $(G4INSTALL)
#checkinc:
#	ls -l $(G4BASE)/geometry/management/include

include $(G4INSTALL)/config/binmake.gmk
