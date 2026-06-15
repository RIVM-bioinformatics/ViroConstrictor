FROM mambaorg/micromamba:latest

COPY ./ViroConstrictor/workflow/envs/mr_scripts.yaml /install.yml

LABEL org.opencontainers.image.description="Supplementary scripts for the ViroConstrictor MR workflow."

USER root

ARG UID=10001

# Create a non-login runtime user without apt-get dependencies.
RUN useradd \
    --home-dir "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser

# Combine micromamba installs and clean in one layer.
RUN micromamba install -q -y -n base git -c conda-forge && \
    micromamba install -q -y -n base -f /install.yml && \
    micromamba clean -q --all --yes

USER appuser

ENV PATH=/opt/conda/bin:$PATH

CMD ["@"]
