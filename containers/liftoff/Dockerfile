# Liftoff 1.6.1. The bioconda package depends on minimap2, which both the
# interval alignment and Liftoff itself call. PATH is set explicitly because
# Nextflow's singularity exec does not run the micromamba entrypoint.
FROM mambaorg/micromamba:1.5.10-noble
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes
ENV PATH="/opt/conda/bin:${PATH}"
ENTRYPOINT []
CMD ["/bin/bash"]
