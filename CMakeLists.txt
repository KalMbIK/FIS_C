cmake_minimum_required(VERSION 3.9)
project(FIS_C C)

set(CMAKE_C_STANDARD 99)

add_compile_options(-Wall -O3)

set(MMio lib/mmio.h lib/mmio.c)
set(wp1 src/suggested_work1.c src/suggested_work1.h)
set(wp2 src/suggested_work2.c src/suggested_work2.h)
set(MyNP legacy/myMatrices.h legacy/myMatrices.c legacy/myNumPyV2.h legacy/myNumPyV2.c)
set(main src/main.c)

add_executable(FIS_C ${MMio} ${MyNP} ${wp1} ${wp2} ${main})
#add_executable(FIS_C ${MMio} main.c)