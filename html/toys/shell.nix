{ pkgs ? import <nixpkgs> {} }:
pkgs.llvmPackages.stdenv.mkDerivation {
  name = "polyhedra dev env";
  stdenv = pkgs.clangStdenv;
  nativeBuildInputs = with pkgs.buildPackages; [
    cmake
    emscripten
    wget
  ];
  shellHook = ''
    export EM_CACHE=${builtins.toString ./.}/.cache
  '';
}
