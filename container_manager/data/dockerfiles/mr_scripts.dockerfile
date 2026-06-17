# Template for the dockerfiles that will be generated.
# Dont change the name of this file, even though it is not very descriptive.
# When the name is Dockerfile, vscode will assume its a dockerfile even with a .j2 extension and provide syntax highlighting and linting.
FROM mambaorg/micromamba:git-d12544c-cuda11.7.1-ubuntu22.04

COPY ./ViroConstrictor/workflow/envs/mr_scripts.yaml /install.yml

LABEL org.opencontainers.image.description="Supplementary scripts for the ViroConstrictor MR workflow."

USER root

ARG UID=10001

RUN useradd \
    --home-dir "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser && \
    micromamba install -q -y -n base -f /install.yml && \
    micromamba clean -q --all --yes

USER appuser

ENV PATH=/opt/conda/bin:$PATH

CMD ["@"]
