# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A

# Include any dependencies generated for this target.
include CMakeFiles/project2A.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/project2A.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/project2A.dir/flags.make

CMakeFiles/project2A.dir/project2A.cxx.o: CMakeFiles/project2A.dir/flags.make
CMakeFiles/project2A.dir/project2A.cxx.o: project2A.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/project2A.dir/project2A.cxx.o"
	/Library/Developer/CommandLineTools/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/project2A.dir/project2A.cxx.o -c /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A/project2A.cxx

CMakeFiles/project2A.dir/project2A.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/project2A.dir/project2A.cxx.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A/project2A.cxx > CMakeFiles/project2A.dir/project2A.cxx.i

CMakeFiles/project2A.dir/project2A.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/project2A.dir/project2A.cxx.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A/project2A.cxx -o CMakeFiles/project2A.dir/project2A.cxx.s

# Object files for target project2A
project2A_OBJECTS = \
"CMakeFiles/project2A.dir/project2A.cxx.o"

# External object files for target project2A
project2A_EXTERNAL_OBJECTS =

project2A.app/Contents/MacOS/project2A: CMakeFiles/project2A.dir/project2A.cxx.o
project2A.app/Contents/MacOS/project2A: CMakeFiles/project2A.dir/build.make
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkDomainsChemistry-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersFlowPaths-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersGeneric-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersHyperTree-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersParallelImaging-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersProgrammable-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersSMP-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersSelection-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersTexture-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersVerdict-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkverdict-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkGeovisCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkproj4-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOAMR-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOEnSight-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOExodus-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOExport-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingGL2PS-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingContextOpenGL-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkgl2ps-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOImport-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOInfovis-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtklibxml2-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOLSDyna-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOMINC-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOMovie-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkoggtheora-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOPLY-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOParallel-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkjsoncpp-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOParallelXML-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOSQL-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtksqlite-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOVideo-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingMath-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingMorphological-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingStatistics-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingStencil-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkInteractionImage-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingImage-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingLIC-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingLOD-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingVolumeOpenGL-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkViewsContext2D-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkViewsInfovis-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersAMR-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersParallel-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkexoIIc-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIONetCDF-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkNetCDF_cxx-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkNetCDF-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkhdf5_hl-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkhdf5-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkParallelCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOXML-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOGeometry-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOXMLParser-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkexpat-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOLegacy-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingOpenGL-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkChartsCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingContext2D-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersImaging-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkInfovisLayout-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkInfovisCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkViewsCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkInteractionWidgets-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersHybrid-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingGeneral-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingSources-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersModeling-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingHybrid-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOImage-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkDICOMParser-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkIOCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkmetaio-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkpng-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtktiff-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkjpeg-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkInteractionStyle-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingAnnotation-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingColor-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingVolume-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingLabel-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingFreeType-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkRenderingCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkCommonColor-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersExtraction-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersStatistics-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingFourier-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkImagingCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkalglib-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersGeometry-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersSources-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersGeneral-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkFiltersCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkCommonExecutionModel-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkCommonComputationalGeometry-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkCommonDataModel-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkCommonMisc-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkCommonTransforms-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkCommonMath-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkCommonSystem-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkCommonCore-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtksys-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkftgl-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkfreetype-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: /Users/alonzoaltamirano/Repos/BUILD-VTK-OLD/lib/libvtkzlib-6.3.1.dylib
project2A.app/Contents/MacOS/project2A: CMakeFiles/project2A.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable project2A.app/Contents/MacOS/project2A"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/project2A.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/project2A.dir/build: project2A.app/Contents/MacOS/project2A

.PHONY : CMakeFiles/project2A.dir/build

CMakeFiles/project2A.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/project2A.dir/cmake_clean.cmake
.PHONY : CMakeFiles/project2A.dir/clean

CMakeFiles/project2A.dir/depend:
	cd /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A /Users/alonzoaltamirano/Documents/UO_Coursework/CIS441/project2A/CMakeFiles/project2A.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/project2A.dir/depend
