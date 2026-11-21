# SINEderella container (Docker). The same image as containers/sinederella.def.
# Build from the repository root:
#   docker build -f containers/Dockerfile --build-arg PIPELINE_COMMIT=$(git describe --always --dirty) -t sinederella .
# Run (the working directory holds the genome and the consensus bank; the run directory is written there):
#   docker run --rm -u $(id -u):$(id -g) -v "$PWD":/work -w /work sinederella SINEderella genome.fa consensuses.fa
#   docker run --rm -u $(id -u):$(id -g) -v "$PWD":/work -w /work sinederella publish_run.sh <RUN_DIR> <SPECIES_CODE>
FROM mambaorg/micromamba:2.3.0-debian12-slim

ARG PIPELINE_COMMIT=unknown

# besides environment.yml, two things a server has in its system and the base image does not:
# gcc: SubFam compiles its k-mer ordering program at run time (the python fallback gives the same order);
# util-linux: column (SINEderella_multi summary table)
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /opt/SINEderella/environment.yml
RUN micromamba install -y -n base -f /opt/SINEderella/environment.yml gcc util-linux \
 && micromamba clean --all --yes \
 && rm -rf /opt/conda/pkgs

# .dockerignore keeps .git and caches out
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/SINEderella
# SINEderella writes this into every run's manifest.txt (there is no git inside the image)
RUN echo "${PIPELINE_COMMIT}" > /opt/SINEderella/PIPELINE_COMMIT \
 && chmod -R a+rX /opt/SINEderella

ENV PATH=/opt/SINEderella:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    LANG=C.UTF-8 \
    PYTHONNOUSERSITE=1 \
    SINEDERELLA_TOOLCHAIN= \
    MPLCONFIGDIR=/tmp/matplotlib

LABEL sinederella.repository="https://github.com/Toki-bio/SINEderella" \
      sinederella.commit="${PIPELINE_COMMIT}"

# MPLCONFIGDIR: with -u <uid> there is no writable home for matplotlib's cache
# no activation script is needed: the tools are found through PATH
ENTRYPOINT []
CMD ["SINEderella"]
