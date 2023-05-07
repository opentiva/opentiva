{
  inputs = {
    # https://github.com/DavHau/mach-nix
    mach-nix.url = "github:DavHau/mach-nix";
    pypi.url = "github:DavHau/pypi-deps-db";
    pypi.flake = false;
    mach-nix.inputs.pypi-deps-db.follows = "pypi";
  };

  outputs = {self, nixpkgs, mach-nix, pypi }@inp:
    let
      l = nixpkgs.lib // builtins;
      supportedSystems = [ "x86_64-linux" "aarch64-darwin" ];
      forAllSystems = f: l.genAttrs supportedSystems
        (system: f system (import nixpkgs {inherit system;}));
    in
    {
      # enter this python environment by executing `nix shell .`
      defaultPackage = forAllSystems (system: pkgs: mach-nix.lib."${system}".mkPython {
        python = "python39Full";
        requirements = builtins.readFile ./requirements.txt + builtins.readFile ./optional_requirements.txt;
      });
    };
}

