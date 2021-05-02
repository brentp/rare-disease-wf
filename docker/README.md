# Dockerfile.manta

This is a simple manta install with ini file adjusted to increase sensitivity.
It also includes [svimmer](https://github.com/DecodeGenetics/svimmer) and [graphtyper](https://github.com/DecodeGenetics/graphtyper) to
merge and then joint-call svs. The [graphtyper sv paper](https://www.nature.com/articles/s41467-019-13341-9) nicely shows that this combination
gives great performance for SVs.

This is currently available on dockerhub as: brentp/manta-graphtyper:v0.0.3

