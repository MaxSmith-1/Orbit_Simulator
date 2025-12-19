Steps for Compiling on WSL:

0. Download mingw-32 / cmake

1. Download conan: ```pip install conan```

2. Use conan to install dependencies: ```conan install .. --output-folder=. --build=missing --profile:build=../mingw_profile --profile:host=../mingw_profile```

2. Compile: ```cmake ..``` then ```make -j``` 

3. Cd into root directory and to run ```./executables/Orbit_Simulator```

Demo Run: ```Orbit_Simulator.exe 10110.447864 config/Spacecraft/demo_spacecraft.json config/Bodies/earth.json```
