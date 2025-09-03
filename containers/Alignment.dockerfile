FROM mambaorg/micromamba:latest

COPY ./ViroConstrictor/workflow/envs/Alignment.yaml /install.yml
COPY ./ViroConstrictor/workflow/files/ /files/
COPY ./ViroConstrictor/workflow/wrappers/ /wrappers/

LABEL org.opencontainers.image.description="Sequence alignment processes and tools for the ViroConstrictor workflow."

USER root
RUN apt-get update && apt-get install -y adduser

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
