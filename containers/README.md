# ViroConstrictor Containers

This directory contains everything needed to define, generate, build, test, and publish workflow containers used by ViroConstrictor.

## What Is In This Folder

- `Dockerfile.j2`: Shared Jinja2 template for all container Dockerfiles.
- `dockerfiles.yaml`: Variables and per-container settings used by the template.
- `generate_dockerfiles.py`: Renders Dockerfiles from the template and YAML config.
- `dockerfiles/*.dockerfile`: Generated Dockerfiles used for building container images.
- `build_containers.py`: Builds missing container images as OCI tar archives.
- `add_OCI_to_docker_engine.py`: Loads built OCI tar files into Docker.
- `tag_and_push_containers.py`: Tags and pushes built images to GHCR.
- `convert_artifact_containers_for_apptainer.py`: Converts OCI tar artifacts to `.sif`.
- `pull_published_containers.py`: Pulls published containers for testing.

## Source Of Truth

Dockerfiles are generated, not hand-maintained.

1. Edit `dockerfiles.yaml` (and optionally `Dockerfile.j2`).
2. Run:

```bash
python containers/generate_dockerfiles.py
```

3. Commit the template/config and regenerated `dockerfiles/*.dockerfile` files.

## How Versioning Works (Hash-Based Tags)

Container version tags are derived from the SHA-256 hash of each workflow environment YAML in:

- `ViroConstrictor/workflow/envs/*.yaml`

The helper function `fetch_hashes()` computes a 6-character hash per env file. That hash becomes part of the image tag.

Examples:

- `viroconstrictor_alignment:<hash>`
- `viroconstrictor_clean:<hash>`

At runtime, Snakemake rules reference `.sif` names with the same hash pattern:

- `viroconstrictor_alignment_<hash>.sif`

## Conda Vs Containers In The Pipeline

Rules define both `conda:` and `container:` directives.

Global execution mode is selected in user config (`REPRODUCTION.repro_method`):

- `containers` -> Snakemake runs with Apptainer/Singularity.
- `conda` -> Snakemake runs with Conda environments.

This switch is configured in `ViroConstrictor/workflow_config.py`.

## Runtime Container Resolution

When `repro_method=containers`:

1. ViroConstrictor determines required containers from current env hashes.
2. It checks local cache (`REPRODUCTION.container_cache_path`) for matching `.sif` files.
3. Missing images are pulled from GHCR using Apptainer/Singularity.

Registry namespace used:

- `ghcr.io/rivm-bioinformatics`

## CI Build, Test, and Publish Flow

### 1) Build and Test Workflow

Defined in `.github/workflows/build_and_test.yml`.

- Runs `python containers/build_containers.py`.
- Skips builds for tags that already exist upstream.
- Writes `containers/builtcontainers.json`.
- Saves built OCI images as `.tar` and uploads them as an artifact.
- Converts built OCI containers to `.sif` for testing.
- Pulls already published containers to fill gaps.

### 2) Publish Workflow

Defined in `.github/workflows/publish_containers.yaml`.

- Triggered on release publication (or manual dispatch).
- Downloads artifact from build workflow.
- Loads OCI tar images into Docker.
- Tags each image to GHCR.
- Pushes images to `ghcr.io/rivm-bioinformatics`.

## Local Build Utility

For local development, `build_local_containers.py` can:

1. Build OCI images.
2. Convert them to `.sif`.
3. Move `.sif` files into your configured `container_cache_path`.

## Notes

- Current image invalidation depends on env YAML hashes, not arbitrary Python script changes.
- `Dockerfile.j2` uses `useradd` for runtime user creation (no `apt-get install adduser`).