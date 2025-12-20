Steps for Compiling on Windows:


0. Download mingw-32 / cmake

1. Download conan (pip install conan)

2. In build folder, run (conan install .. --output-folder=. --build=missing --profile:build=../mingw_profile --profile:host=../mingw_profile
)

2. Still in build folder, run (cmake ..) then (mingw32-make) 

3. Cd into executable folder and run (Orbit_Simulator.exe)

Demo Run: Orbit_Simulator.exe 20110.447864 config/Spacecraft/demo_spacecraft.json config/Bodies/earth.json 