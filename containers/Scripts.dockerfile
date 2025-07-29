FROM mambaorg/micromamba:latest

COPY ./ViroConstrictor/workflow/envs/Scripts.yaml /install.yml
COPY ./ViroConstrictor/workflow/files/ /files/
COPY ./ViroConstrictor/workflow/wrappers/ /wrappers/
COPY ./ViroConstrictor/workflow/scripts/ /scripts/

LABEL org.opencontainers.image.description="Supplementary scripts for the ViroConstrictor workflow."

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

RUN touch /__init__.py
    
RUN micromamba install -q -y -n base git -c conda-forge && \
    micromamba install -q -y -n base -f /install.yml && \
    micromamba clean -q --all --yes

USER appuser

ENV PATH=/opt/conda/bin:$PATH

CMD ["@"]
