# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

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
CMAKE_COMMAND = /opt/local/bin/cmake

# The command to remove a file.
RM = /opt/local/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/angus/CPP_testing/firstCode/src

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/angus/CPP_testing/firstCode/build

# Include any dependencies generated for this target.
include CMakeFiles/tdb.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/tdb.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tdb.dir/flags.make

CMakeFiles/tdb.dir/main.cpp.o: CMakeFiles/tdb.dir/flags.make
CMakeFiles/tdb.dir/main.cpp.o: /Users/angus/CPP_testing/firstCode/src/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/angus/CPP_testing/firstCode/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tdb.dir/main.cpp.o"
	/opt/local/bin/x86_64-apple-darwin13-g++-mp-4.9   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tdb.dir/main.cpp.o -c /Users/angus/CPP_testing/firstCode/src/main.cpp

CMakeFiles/tdb.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tdb.dir/main.cpp.i"
	/opt/local/bin/x86_64-apple-darwin13-g++-mp-4.9  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/angus/CPP_testing/firstCode/src/main.cpp > CMakeFiles/tdb.dir/main.cpp.i

CMakeFiles/tdb.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tdb.dir/main.cpp.s"
	/opt/local/bin/x86_64-apple-darwin13-g++-mp-4.9  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/angus/CPP_testing/firstCode/src/main.cpp -o CMakeFiles/tdb.dir/main.cpp.s

CMakeFiles/tdb.dir/main.cpp.o.requires:

.PHONY : CMakeFiles/tdb.dir/main.cpp.o.requires

CMakeFiles/tdb.dir/main.cpp.o.provides: CMakeFiles/tdb.dir/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/tdb.dir/build.make CMakeFiles/tdb.dir/main.cpp.o.provides.build
.PHONY : CMakeFiles/tdb.dir/main.cpp.o.provides

CMakeFiles/tdb.dir/main.cpp.o.provides.build: CMakeFiles/tdb.dir/main.cpp.o


CMakeFiles/tdb.dir/jsoncpp.cpp.o: CMakeFiles/tdb.dir/flags.make
CMakeFiles/tdb.dir/jsoncpp.cpp.o: /Users/angus/CPP_testing/firstCode/src/jsoncpp.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/angus/CPP_testing/firstCode/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/tdb.dir/jsoncpp.cpp.o"
	/opt/local/bin/x86_64-apple-darwin13-g++-mp-4.9   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/tdb.dir/jsoncpp.cpp.o -c /Users/angus/CPP_testing/firstCode/src/jsoncpp.cpp

CMakeFiles/tdb.dir/jsoncpp.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tdb.dir/jsoncpp.cpp.i"
	/opt/local/bin/x86_64-apple-darwin13-g++-mp-4.9  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/angus/CPP_testing/firstCode/src/jsoncpp.cpp > CMakeFiles/tdb.dir/jsoncpp.cpp.i

CMakeFiles/tdb.dir/jsoncpp.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tdb.dir/jsoncpp.cpp.s"
	/opt/local/bin/x86_64-apple-darwin13-g++-mp-4.9  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/angus/CPP_testing/firstCode/src/jsoncpp.cpp -o CMakeFiles/tdb.dir/jsoncpp.cpp.s

CMakeFiles/tdb.dir/jsoncpp.cpp.o.requires:

.PHONY : CMakeFiles/tdb.dir/jsoncpp.cpp.o.requires

CMakeFiles/tdb.dir/jsoncpp.cpp.o.provides: CMakeFiles/tdb.dir/jsoncpp.cpp.o.requires
	$(MAKE) -f CMakeFiles/tdb.dir/build.make CMakeFiles/tdb.dir/jsoncpp.cpp.o.provides.build
.PHONY : CMakeFiles/tdb.dir/jsoncpp.cpp.o.provides

CMakeFiles/tdb.dir/jsoncpp.cpp.o.provides.build: CMakeFiles/tdb.dir/jsoncpp.cpp.o


# Object files for target tdb
tdb_OBJECTS = \
"CMakeFiles/tdb.dir/main.cpp.o" \
"CMakeFiles/tdb.dir/jsoncpp.cpp.o"

# External object files for target tdb
tdb_EXTERNAL_OBJECTS =

tdb: CMakeFiles/tdb.dir/main.cpp.o
tdb: CMakeFiles/tdb.dir/jsoncpp.cpp.o
tdb: CMakeFiles/tdb.dir/build.make
tdb: CMakeFiles/tdb.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/angus/CPP_testing/firstCode/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable tdb"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tdb.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tdb.dir/build: tdb

.PHONY : CMakeFiles/tdb.dir/build

CMakeFiles/tdb.dir/requires: CMakeFiles/tdb.dir/main.cpp.o.requires
CMakeFiles/tdb.dir/requires: CMakeFiles/tdb.dir/jsoncpp.cpp.o.requires

.PHONY : CMakeFiles/tdb.dir/requires

CMakeFiles/tdb.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tdb.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tdb.dir/clean

CMakeFiles/tdb.dir/depend:
	cd /Users/angus/CPP_testing/firstCode/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/angus/CPP_testing/firstCode/src /Users/angus/CPP_testing/firstCode/src /Users/angus/CPP_testing/firstCode/build /Users/angus/CPP_testing/firstCode/build /Users/angus/CPP_testing/firstCode/build/CMakeFiles/tdb.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tdb.dir/depend
