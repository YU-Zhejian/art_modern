# Readme for Packing

The automated packing system will do the following:

- Build Singularity containers for builders.
- Execute the builder using the source code.
- Package the output artifacts.

## Packing Issues due to Dependencies

- AlmaLinux 8: The HTSLib shipped is 1.9, lower than required 1.14.
- AlmaLinux 9: lacking HTSLib.
- Ubuntu 22.04: The HTSLib shipped is 1.13, lower than required 1.14.

Users of those distributions should install required libraries manually or use bundled versions and build from source.
