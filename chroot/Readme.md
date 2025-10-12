# Readme for `chroot`

Here contains some `chroot`s for testing or building purposes. Password-less access to root privilidge via (sudo) is required for most of the operations.

Use `make all` to automatically build all `chroot` images.

`chroot`s:

- Ubuntu 18.04 LTS.
- Debian Oldoldstable.

**NOTE** We would currently not use Ubuntu 18.04 LTS and Debian Oldoldstable for building official DEB packages as they require the inclusion of `BS::thread_pool` which does not comply with the policy of Debian and Ubuntu.
