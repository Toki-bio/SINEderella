# SINEderella container (Docker). The same image as containers/sinederella.def.
# Build from the repository root:
#   docker build -f containers/Dockerfile --build-arg PIPELINE_COMMIT=$(git describe --always --dirty) -t sinederella .
# Run (the working directory holds the genome and the consensus bank; the run directory is written there):
#   docker run --rm -u $(id -u):$(id -g) -v "$PWD":/work -w /work sinederella SINEderella genome.fa consensuses.fa
#   docker run --rm -u $(id -u):$(id -g) -v "$PWD":/work -w /work sinederella publish_run.sh <RUN_DIR> <SPECIES_CODE>
FROM mambaorg/micromamba:2.3.0-debian12-slim

ARG PIPELINE_COMMIT=unknown

USER root
# gcc: SubFam compiles its k-mer ordering program at run time (python fallback gives the same order);
# bsdextrautils: column (SINEderella_multi summary table)
RUN apt-get update \
 && apt-get install -y --no-install-recommends gcc libc6-dev bsdextrautils \
 && rm -rf /var/lib/apt/lists/*

COPY environment.yml /opt/SINEderella/environment.yml
RUN micromamba install -y -n base -f /opt/SINEderella/environment.yml \
 && micromamba clean --all --yes

# .dockerignore keeps .git, caches and run directories out
COPY . /opt/SINEderella
# SINEderella writes this into every run's manifest.txt (there is no git inside the image)
RUN echo "${PIPELINE_COMMIT}" > /opt/SINEderella/PIPELINE_COMMIT \
 && chmod -R a+rX /opt/SINEderella /opt/conda

ENV PATH=/opt/SINEderella:/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    LANG=C.UTF-8 \
    PYTHONNOUSERSITE=1 \
    SINEDERELLA_TOOLCHAIN=

LABEL org.opencontainers.image.source="https://github.com/Toki-bio/SINEderella" \
      org.opencontainers.image.title="SINEderella"

# no activation script is needed: the tools are found through PATH
ENTRYPOINT []
CMD ["SINEderella"]
