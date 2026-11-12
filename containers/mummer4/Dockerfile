# MUMmer 4.0.1 (nucmer, mummerplot) and gnuplot 5.4.
# Nextflow runs singularity exec, which does not invoke the micromamba entrypoint,
# so the env bin directory has to be on PATH by itself.
FROM mambaorg/micromamba:1.5.10-noble
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes
ENV PATH="/opt/conda/bin:${PATH}"
ENTRYPOINT []
CMD ["/bin/bash"]
