#!/usr/bin/sh
case $1 in
"load") ;;
"clean")
    rm -rf build
    mkdir build
    rm -rf bin
    mkdir bin
    ;;
"build")
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER="gcc-11" -DCMAKE_CXX_COMPILER="g++-11" ..
    make -j8
    # add share lib for mpi version
    cd ..
    deplist=$(ldd bin/SiprosEnsembleMPI | awk '{if (match($3,"/")){ print $3}}')
    mkdir bin/libSiprosEnsembleMPI
    cp -L -n $deplist bin/libSiprosEnsembleMPI
    ;;
"buildTick")
    cd build
    cmake -DCMAKE_BUILD_TYPE=Release -DTicktock=Ticktock
    make -j4
    ;;
"debug")
    cd build
    cmake -DCMAKE_BUILD_TYPE=Debug
    make -j4
    ;;
"make")
    cd build
    make
    ;;
"run")
    cd timeCompare
    starttime=$(date +'%Y-%m-%d %H:%M:%S')
    ../bin/SiprosV3omp -c ../../data/SiproConfig.N15_0Pct.cfg \
        -f ../../data/AMD_DynamicSIP_SampleD_TimePoint0_BRmixed_WC_Velos_OrbiMS2_Run2_020210_09.FT2 \
        -o . -s
    endtime=$(date +'%Y-%m-%d %H:%M:%S')
    start_seconds=$(date --date="$starttime" +%s)
    end_seconds=$(date --date="$endtime" +%s)
    echo "running time： "$((end_seconds - start_seconds))"s"
    ;;
*)
    ./make "build"
    ;;
esac
