# NanoLooper

## Instructions

Log in to UAF machine and create a work directory

    ssh uaf-10.t2.ucsd.edu
    git clone --recursive https://github.com/sgnoohc/NanoLooper.git
    cd NanoLooper

Set up the environment

    source setup.sh

Compile the code

    make -j # Compile

Run the looper

    sh run.sh
    sh hadd.sh
    python plot.py

Examine the code

    vim NanoLooper/main.cc


