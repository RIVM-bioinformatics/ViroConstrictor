FROM mambaorg/micromamba:latest

COPY ./ViroConstrictor/workflow/envs/ORF_analysis.yaml /install.yml

LABEL org.opencontainers.image.description="ORF analysis processes and tools for the ViroConstrictor workflow."

USER root

ARG UID=10001

RUN useradd --no-create-home --shell /sbin/nologin --uid "${UID}" appuser && \
    micromamba install -q -y -n base -f /install.yml && \
    micromamba clean -q --all --yes && \
    rm -rf /opt/conda/pkgs && \
    rm -rf /opt/conda/x86_64-conda-linux-gnu && \
    rm -rf /opt/conda/include && \
    rm -rf /opt/conda/share/doc /opt/conda/share/man /opt/conda/share/info

USER appuser

ENV PATH=/opt/conda/bin:$PATH

CMD ["@"]
