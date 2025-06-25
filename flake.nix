{
  description = "VTK and C++ development setup";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-24.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};
      in
      {
        devShells.default = pkgs.mkShell {
          buildInputs =
            with pkgs;
            [
              cmake
              pkg-config
              gcc

              vtk
              libGL
              glew

              gdb
            ];

          shellHook = ''
            export LSAN_OPTIONS=suppressions=$PWD/suppr.txt
          '';
        };
      }
    );
}
