FROM mambaorg/micromamba:latest

COPY ./ViroConstrictor/workflow/envs/ORF_analysis.yaml /install.yml

LABEL org.opencontainers.image.description="ORF analysis processes and tools for the ViroConstrictor workflow."

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
