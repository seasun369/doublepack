# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.21

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /home/anaconda3/lib/python3.9/site-packages/cmake/data/bin/cmake

# The command to remove a file.
RM = /home/anaconda3/lib/python3.9/site-packages/cmake/data/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/extra_space/ztx/turbopack/secure-computation-library

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release

# Include any dependencies generated for this target.
include CMakeFiles/scl.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/scl.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/scl.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/scl.dir/flags.make

CMakeFiles/scl.dir/src/scl/prg.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/prg.cc.o: ../../src/scl/prg.cc
CMakeFiles/scl.dir/src/scl/prg.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/scl.dir/src/scl/prg.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/prg.cc.o -MF CMakeFiles/scl.dir/src/scl/prg.cc.o.d -o CMakeFiles/scl.dir/src/scl/prg.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/prg.cc

CMakeFiles/scl.dir/src/scl/prg.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/prg.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/prg.cc > CMakeFiles/scl.dir/src/scl/prg.cc.i

CMakeFiles/scl.dir/src/scl/prg.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/prg.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/prg.cc -o CMakeFiles/scl.dir/src/scl/prg.cc.s

CMakeFiles/scl.dir/src/scl/hash.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/hash.cc.o: ../../src/scl/hash.cc
CMakeFiles/scl.dir/src/scl/hash.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/scl.dir/src/scl/hash.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/hash.cc.o -MF CMakeFiles/scl.dir/src/scl/hash.cc.o.d -o CMakeFiles/scl.dir/src/scl/hash.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/hash.cc

CMakeFiles/scl.dir/src/scl/hash.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/hash.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/hash.cc > CMakeFiles/scl.dir/src/scl/hash.cc.i

CMakeFiles/scl.dir/src/scl/hash.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/hash.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/hash.cc -o CMakeFiles/scl.dir/src/scl/hash.cc.s

CMakeFiles/scl.dir/src/scl/math/str.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/math/str.cc.o: ../../src/scl/math/str.cc
CMakeFiles/scl.dir/src/scl/math/str.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/scl.dir/src/scl/math/str.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/math/str.cc.o -MF CMakeFiles/scl.dir/src/scl/math/str.cc.o.d -o CMakeFiles/scl.dir/src/scl/math/str.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/math/str.cc

CMakeFiles/scl.dir/src/scl/math/str.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/math/str.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/math/str.cc > CMakeFiles/scl.dir/src/scl/math/str.cc.i

CMakeFiles/scl.dir/src/scl/math/str.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/math/str.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/math/str.cc -o CMakeFiles/scl.dir/src/scl/math/str.cc.s

CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.o: ../../src/scl/math/fields/mersenne61.cc
CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.o -MF CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.o.d -o CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/math/fields/mersenne61.cc

CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/math/fields/mersenne61.cc > CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.i

CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/math/fields/mersenne61.cc -o CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.s

CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.o: ../../src/scl/math/fields/mersenne127.cc
CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.o -MF CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.o.d -o CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/math/fields/mersenne127.cc

CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/math/fields/mersenne127.cc > CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.i

CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/math/fields/mersenne127.cc -o CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.s

CMakeFiles/scl.dir/src/scl/net/config.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/net/config.cc.o: ../../src/scl/net/config.cc
CMakeFiles/scl.dir/src/scl/net/config.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/scl.dir/src/scl/net/config.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/net/config.cc.o -MF CMakeFiles/scl.dir/src/scl/net/config.cc.o.d -o CMakeFiles/scl.dir/src/scl/net/config.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/config.cc

CMakeFiles/scl.dir/src/scl/net/config.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/net/config.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/config.cc > CMakeFiles/scl.dir/src/scl/net/config.cc.i

CMakeFiles/scl.dir/src/scl/net/config.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/net/config.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/config.cc -o CMakeFiles/scl.dir/src/scl/net/config.cc.s

CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.o: ../../src/scl/net/mem_channel.cc
CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.o -MF CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.o.d -o CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/mem_channel.cc

CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/mem_channel.cc > CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.i

CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/mem_channel.cc -o CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.s

CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.o: ../../src/scl/net/tcp_channel.cc
CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.o -MF CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.o.d -o CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/tcp_channel.cc

CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/tcp_channel.cc > CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.i

CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/tcp_channel.cc -o CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.s

CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.o: ../../src/scl/net/threaded_sender.cc
CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.o -MF CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.o.d -o CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/threaded_sender.cc

CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/threaded_sender.cc > CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.i

CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/threaded_sender.cc -o CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.s

CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.o: ../../src/scl/net/tcp_utils.cc
CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.o -MF CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.o.d -o CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/tcp_utils.cc

CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/tcp_utils.cc > CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.i

CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/tcp_utils.cc -o CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.s

CMakeFiles/scl.dir/src/scl/net/network.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/net/network.cc.o: ../../src/scl/net/network.cc
CMakeFiles/scl.dir/src/scl/net/network.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/scl.dir/src/scl/net/network.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/net/network.cc.o -MF CMakeFiles/scl.dir/src/scl/net/network.cc.o.d -o CMakeFiles/scl.dir/src/scl/net/network.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/network.cc

CMakeFiles/scl.dir/src/scl/net/network.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/net/network.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/network.cc > CMakeFiles/scl.dir/src/scl/net/network.cc.i

CMakeFiles/scl.dir/src/scl/net/network.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/net/network.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/network.cc -o CMakeFiles/scl.dir/src/scl/net/network.cc.s

CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.o: ../../src/scl/net/discovery/server.cc
CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.o -MF CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.o.d -o CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/discovery/server.cc

CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/discovery/server.cc > CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.i

CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/discovery/server.cc -o CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.s

CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.o: CMakeFiles/scl.dir/flags.make
CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.o: ../../src/scl/net/discovery/client.cc
CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.o: CMakeFiles/scl.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.o -MF CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.o.d -o CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.o -c /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/discovery/client.cc

CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/discovery/client.cc > CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.i

CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/extra_space/ztx/turbopack/secure-computation-library/src/scl/net/discovery/client.cc -o CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.s

# Object files for target scl
scl_OBJECTS = \
"CMakeFiles/scl.dir/src/scl/prg.cc.o" \
"CMakeFiles/scl.dir/src/scl/hash.cc.o" \
"CMakeFiles/scl.dir/src/scl/math/str.cc.o" \
"CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.o" \
"CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.o" \
"CMakeFiles/scl.dir/src/scl/net/config.cc.o" \
"CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.o" \
"CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.o" \
"CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.o" \
"CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.o" \
"CMakeFiles/scl.dir/src/scl/net/network.cc.o" \
"CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.o" \
"CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.o"

# External object files for target scl
scl_EXTERNAL_OBJECTS =

libscl.so.0.3: CMakeFiles/scl.dir/src/scl/prg.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/hash.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/math/str.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/math/fields/mersenne61.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/math/fields/mersenne127.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/net/config.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/net/mem_channel.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/net/tcp_channel.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/net/threaded_sender.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/net/tcp_utils.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/net/network.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/net/discovery/server.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/src/scl/net/discovery/client.cc.o
libscl.so.0.3: CMakeFiles/scl.dir/build.make
libscl.so.0.3: CMakeFiles/scl.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Linking CXX shared library libscl.so"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/scl.dir/link.txt --verbose=$(VERBOSE)
	$(CMAKE_COMMAND) -E cmake_symlink_library libscl.so.0.3 libscl.so.0.3 libscl.so

libscl.so: libscl.so.0.3
	@$(CMAKE_COMMAND) -E touch_nocreate libscl.so

# Rule to build all files generated by this target.
CMakeFiles/scl.dir/build: libscl.so
.PHONY : CMakeFiles/scl.dir/build

CMakeFiles/scl.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/scl.dir/cmake_clean.cmake
.PHONY : CMakeFiles/scl.dir/clean

CMakeFiles/scl.dir/depend:
	cd /mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/extra_space/ztx/turbopack/secure-computation-library /mnt/extra_space/ztx/turbopack/secure-computation-library /mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release /mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release /mnt/extra_space/ztx/turbopack/secure-computation-library/build/Release/CMakeFiles/scl.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/scl.dir/depend

