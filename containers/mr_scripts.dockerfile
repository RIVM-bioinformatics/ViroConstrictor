FROM mambaorg/micromamba:latest

COPY ./ViroConstrictor/workflow/envs/mr_scripts.yaml /install.yml

LABEL org.opencontainers.image.description="Supplementary scripts for the ViroConstrictor MR workflow."

USER root

ARG UID=10001

# Combine apt-get commands and clean up in the same layer to reduce size
RUN apt-get update && \
    apt-get install -y --no-install-recommends adduser && \
    adduser \
    --disabled-password \
    --gecos "" \
    --home "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Combine micromamba installs and clean in one layer
# Added --prune to remove unused packages
RUN micromamba install -q -y -n base git -c conda-forge && \
    micromamba install -q -y -n base -f /install.yml && \
    micromamba clean -q --all --yes
USER appuser

ENV PATH=/opt/conda/bin:$PATH

CMD ["@"]
