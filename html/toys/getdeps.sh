#/usr/bin/env bash

root_dir=$(pwd)
bootstrap_flag='false'
patches=()

print_usage() {
  echo "
  -b: Download a bootstrap project
  -d: Specify directory to use, default is pwd
  -p: Apply a patch file
  "
}

while getopts 'bd:p:' flag; do
  case ${flag} in
    b) bootstrap_flag='true' ;;
    d) root_dir=${OPTARG} ;;
    p) patches+=( "${OPTARG}" ) ;;
    *) print_usage
       exit 0 ;;
  esac
done

absolute_patches=()
for p in ${patches[@]}; do
  absolute_patches+=($(realpath $p))
done

mkdir -p $root_dir && pushd $root_dir
root_dir=$(pwd)
mkdir -p system

# download magnum and a bootstrap project

if [ "$bootstrap_flag" = true ]
then
  wget https://github.com/mosra/magnum-bootstrap/archive/refs/heads/base-emscripten.zip -O base.zip
  bsdtar --strip-components=1 -xf base.zip -C .
  rm base.zip
fi

wget https://github.com/mosra/toolchains/archive/refs/heads/master.zip -O toolchains.zip
mkdir -p toolchains
bsdtar --strip-components=1 -xf toolchains.zip -C toolchains
rm toolchains.zip

wget https://github.com/mosra/corrade/archive/refs/heads/master.zip -O corrade.zip
mkdir -p corrade
bsdtar --strip-components=1 -xf corrade.zip -C corrade
rm corrade.zip

wget https://github.com/mosra/magnum/archive/refs/heads/master.zip -O magnum.zip
mkdir -p magnum
bsdtar --strip-components=1 -xf magnum.zip -C magnum
rm magnum.zip

# wget https://github.com/libsdl-org/SDL/releases/download/release-2.0.10/SDL2-2.0.10.zip -O SDL2.zip
# mkdir -p SDL2
# bsdtar --strip-components=1 -xf SDL2.zip -C SDL2
# rm SDL2.zip

# apply patches

for p in ${absolute_patches[@]}; do
  patch -p0 < $p
done

# build corrade

pushd corrade

mkdir build && pushd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="${root_dir}/system" \
  -DCMAKE_INSTALL_PREFIX="${root_dir}/system"
cmake --build .
popd

mkdir build-emscripten-wasm && pushd build-emscripten-wasm
cmake .. \
  -DCMAKE_TOOLCHAIN_FILE="${root_dir}/toolchains/generic/Emscripten-wasm.cmake" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="${root_dir}/system" \
  -DCMAKE_INSTALL_PREFIX="${root_dir}/system" \
  -DCORRADE_RC_EXECUTABLE="../build/bin/corrande-rc"
cmake --build .
cmake --build . --target install
popd

popd

# build magnum

pushd magnum

mkdir build-emscripten-wasm && pushd build-emscripten-wasm
cmake .. \
  -DCMAKE_TOOLCHAIN_FILE="${root_dir}/toolchains/generic/Emscripten-wasm.cmake" \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="${root_dir}/system" \
  -DCMAKE_INSTALL_PREFIX="${root_dir}/system" \
  -DCORRADE_RC_EXECUTABLE="${root_dir}/corrade/build/Release/bin/corrade-rc" \
  -DMAGNUM_WITH_EMSCRIPTENAPPLICATION=ON
cmake --build .
cmake --build . --target install
popd

popd

popd
