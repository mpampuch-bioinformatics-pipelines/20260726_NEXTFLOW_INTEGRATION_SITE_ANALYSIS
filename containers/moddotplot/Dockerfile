# modDotPlot 1.0.0 from PyPI. Bioconda 0.9.9 has no --pairs or --sequence,
# which the dot-plot workflow uses to avoid an all-versus-all of every
# reference chromosome. PATH is set explicitly because Nextflow's singularity
# exec does not run the micromamba entrypoint.
FROM mambaorg/micromamba:1.5.10-noble
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml \
    && micromamba clean --all --yes
ENV PATH="/opt/conda/bin:${PATH}"
ENV MPLCONFIGDIR=/tmp/matplotlib
ENTRYPOINT []
CMD ["/bin/bash"]
