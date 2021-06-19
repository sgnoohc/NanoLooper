# Instruction

Log in to UAF machine and create a work directory

    ssh uaf-10.t2.ucsd.edu
    mkdir nano_analysis
    cd nano_analysis

Downloading some packages

    git clone https://github.com/cmstas/NanoTools.git
    git clone https://github.com/sgnoohc/rooutil.git

Set up the environment

    source rooutil/thisrooutil.sh
    source rooutil/root.sh

Compile the code

    cd NanoTools/NanoCORE
    make -j # Compile


