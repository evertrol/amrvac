- Add a macro flag to the compiler flags: `-DGENDIM

  Best to do this in the relevant arch/*defs file; adding it in the
  project makefile won't build the library with the flag, since the
  library has its own makefile, which itself (again) reads the
  relevant *defs file.

  This macro is already added to the default and debug builds.

- Use `setup.pl -d=3` as usual, and run `make`.

- set `autoconvert = F` in the parameter file. This part has not been
  modified for generic n-D use.

- Set the `block_nx?` to 0 for the dimensions that will not be used;
  thus, starting with `block_nx3`.
