# container_manager

Utilities for generating, building, converting, publishing, and syncing workflow containers.

## Overview

The container manager works as a staged pipeline:

1. `generate`: render Dockerfiles from `config_dockerfiles.yaml` + `docker_template_specific.j2`
2. `build`: build Docker images and export `.tar` artifacts, then write a manifest
3. `convert`: convert `.tar` artifacts to `.sif` using Apptainer
4. `publish`: tag and push built images listed in the manifest
5. `sync-cache`: compute required hash-tagged images and pull missing containers to a local cache
6. `merge-manifests`: merge per-runner manifests into one combined manifest

The `local` command runs a convenience flow: `generate -> build -> convert -> cache sync`.

## Entry points

Run as module:

```bash
python -m container_manager --help
```

Run by package path:

```bash
python container_manager/ --help
```

## Required config

Configuration lives in `container_manager/config_dockerfiles.yaml`.

Required top-level keys:

- `registry`: image registry namespace used to form remote image refs
- `container_package_prefix`: image name prefix used for tag/cache naming
- `oci_labels`: mapping of OCI labels injected into built images
- `dockerfiles`: list of Dockerfile build specs

Each `dockerfiles` entry must contain:

- `output`: output Dockerfile path ending in `.dockerfile`
- `env_file`: conda env file used in the template
- `description`: human-readable image description
- `extra_copies`: list of `{src, dest}` pairs copied into image context

## Common commands

Generate Dockerfiles:

```bash
python -m container_manager generate
```

Plan builds without executing:

```bash
python -m container_manager build --dry-run
```

Build using Dockerfiles from one directory and write artifacts to another:

```bash
python -m container_manager build \
	--dockerfiles-dir ./container_manager/data/dockerfiles \
	--output-dir ./container_manager/data/containers
```

Build one image and write/update manifest:

```bash
python -m container_manager build --only Alignment
```

Convert built tar artifacts to SIF:

```bash
python -m container_manager convert
```

Publish manifest images (planned mode):

```bash
python -m container_manager publish --dry-run
```

Merge matrix build manifests:

```bash
python -m container_manager merge-manifests --input-glob "./artifacts/container-build-manifest-*.json"
```

Sync local cache from registry (planned mode):

```bash
python -m container_manager sync-cache --cache-dir ~/.viroconstrictor/containers --dry-run
```

Run the local end-to-end flow:

```bash
python -m container_manager local \
	--cache-dir ~/.viroconstrictor/containers \
	--dockerfiles-dir ./container_manager/data/dockerfiles \
	--output-dir ./container_manager/data/containers
```

## Notes

- `build` expects Dockerfiles to exist. Run `generate` first unless you already generated Dockerfiles.
- `build` no longer takes a `--template` flag; template hashing uses the default template path from configuration constants.
- `build --dry-run` and `convert --dry-run` still write planned status updates to the manifest.
- `publish --dry-run` does not fail on unbuilt items; these are marked as skipped in manifest status.
- Runtime workflow execution uses Apptainer/Singularity `.sif` images, not Docker images directly.
