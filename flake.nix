{
  description = "MOTCO project development environment flake";
  
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    systems.url = "github:nix-systems/default";
  };
  
  outputs =
    {
      self,
      nixpkgs,
      systems,
    }:
    let
      forEachSystem = nixpkgs.lib.genAttrs (import systems);
    in
    {
      devShells = forEachSystem (
        system:
        let
          pkgs = nixpkgs.legacyPackages.${system};
        in
        {
          default = pkgs.mkShell {
            packages = with pkgs; [
              python311
              uv
              R
              rPackages.InterSIM
              rPackages.mixOmics
            ];

            LD_LIBRARY_PATH = pkgs.lib.makeLibraryPath [ pkgs.stdenv.cc.cc.lib ];

            shellHook = ''
              export UV_PYTHON="${pkgs.python311}/bin/python"
              export UV_PYTHON_DOWNLOADS=never

              uv sync --locked --extra test --extra docs
            '';
          };
        }
      );
    };
}
