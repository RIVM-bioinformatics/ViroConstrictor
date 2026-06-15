FROM mambaorg/micromamba:latest

COPY ./ViroConstrictor/workflow/envs/Clean.yaml /install.yml
COPY ./ViroConstrictor/workflow/main/configs/ /configs/

LABEL org.opencontainers.image.description="Data cleaning processes and tools for the ViroConstrictor workflow."

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
