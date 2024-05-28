. intial_cmake
# . setenv.sh
echo $prefix
rm -rf other_build
rm -rf lock-baryon-pi-psel-32Dfine-*
mkdir other_build
cd other_build
cmake ..
make -j12
cd ..

# for n in {1..10};
for n in {1..30};
do
    echo $n
    yhbatch o_run.sh
    # yhbatch run_4.sh
    # yhbatch run_8.sh
    # yhbatch run_16.sh
    # yhbatch run_32.sh
    # yhbatch run_64.sh
done