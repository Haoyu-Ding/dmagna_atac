FROM mambaorg/micromamba:1.5.8

USER root
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && rm -f /tmp/environment.yml

COPY --chown=$MAMBA_USER:$MAMBA_USER bin /opt/dmagna-atac/bin
RUN chmod +x /opt/dmagna-atac/bin/validate_samplesheet.py

ENV PATH=/opt/conda/bin:$PATH
USER $MAMBA_USER
WORKDIR /work
