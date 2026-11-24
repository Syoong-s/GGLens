# Publishing releases

## Preview the archive

After committing the release changes, list the files Git will export from that commit.

```sh
git archive --format=tar.gz --prefix=FQ_GGLens/ HEAD | tar -tzf -
```

The archive contains `LICENSE`, `code/Makefile`, the headers in `code/include/`, the C++ sources in `code/src/`, the Python modules listed in `FIT_MODULES`, and the C++ runtime selection script. Build from the extracted archive root with `make -C code`. In a development checkout, GNU Make reads `code/GNUmakefile` first. It includes the build rules and adds the targets that use `code/tests/`. The release allowlist exports only the build Makefile. Use the tagged repository to access the tests, documentation, and configuration examples. [Archive allowlist](../.gitattributes)

## Publish a version

After pushing the release commit to GitHub, tag that commit with a new version of the form `vMAJOR.MINOR.PATCH` and push the tag.

```sh
git tag -a v0.1.0 -m "Release v0.1.0"
git push origin v0.1.0
```

Before uploading, the workflow extracts the archive, checks its Makefile and the paths reserved for development files, then builds and cleans the extracted code. It uses the repository's Actions token to create a GitHub Release for the tag and attach `FQ_GGLens-v0.1.0.tar.gz` in this example. Check the run in the Actions tab and the archive in the Releases tab. [Release workflow](../.github/workflows/release.yml)
