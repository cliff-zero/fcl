. intial_cmake
echo $prefix
# rm -rf lock*
# rm -rf ./out/out_*
# rm -rf ./err/err_*
# rm -rf ./result/*
# rm -rf core.*
# rm -rf build
# mkdir build
# cd build
# cmake ..
# make -j12
# cd ..

# for n in {1..1};
for n in {1..4};
do
#     echo $n
#     # yhbatch run.sh
#     # yhbatch run_4.sh
#     # yhbatch run_8.sh
#     # yhbatch run_16.sh
    yhbatch run_32.sh
#     # yhbatch run_64.sh
done