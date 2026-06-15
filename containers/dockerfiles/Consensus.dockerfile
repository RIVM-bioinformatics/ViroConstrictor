FROM mambaorg/micromamba:latest

COPY ./ViroConstrictor/workflow/envs/Consensus.yaml /install.yml

LABEL org.opencontainers.image.description="Consensus sequence generation processes and tools for the ViroConstrictor workflow."

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
