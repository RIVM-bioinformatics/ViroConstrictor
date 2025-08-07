FROM mambaorg/micromamba:latest

COPY ./ViroConstrictor/workflow/envs/Clean.yaml /install.yml
COPY ./ViroConstrictor/workflow/main/configs/ /configs/

LABEL org.opencontainers.image.description="Data cleaning processes and tools for the ViroConstrictor workflow."

USER root

ARG UID=10001
RUN adduser \
    --disabled-password \
    --gecos "" \
    --home "/nonexistent" \
    --shell "/sbin/nologin" \
    --no-create-home \
    --uid "${UID}" \
    appuser

RUN micromamba install -q -y -n base git -c conda-forge && \
    micromamba install -q -y -n base -f /install.yml && \
    micromamba clean -q --all --yes

USER appuser

ENV PATH=/opt/conda/bin:$PATH

CMD ["@"]
