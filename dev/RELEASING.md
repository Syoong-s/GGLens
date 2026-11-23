# Publishing releases

## Preview the archive

After committing the release changes, list the files Git will export from that commit.

```sh
git archive --format=tar.gz --prefix=FQ_GGLens/ HEAD | tar -tzf -
```

The archive contains `LICENSE`, `code/Makefile`, the headers in `code/include/`, the C++ sources in `code/src/`, and the five Python modules listed in `FIT_MODULES`. Build from the extracted archive root with `make -C code`. Update the allowlist when a new build input falls outside its current rules. Use the tagged repository to access the documentation and configuration examples. [Archive allowlist](../.gitattributes)

## Publish a version

After pushing the release commit to GitHub, tag that commit with a new version of the form `vMAJOR.MINOR.PATCH` and push the tag.

```sh
git tag -a v0.1.0 -m "Release v0.1.0"
git push origin v0.1.0
```

The workflow uses the repository's Actions token to create a GitHub Release for the tag and attach `FQ_GGLens-v0.1.0.tar.gz` in this example. Check the run in the Actions tab and the archive in the Releases tab. [Release workflow](../.github/workflows/release.yml)
