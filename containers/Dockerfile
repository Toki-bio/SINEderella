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

# SINEplot (github.com/Toki-bio/SINEplot, fixed commit): the bit-score layout panel of the report; step6 finds it on PATH
RUN /opt/conda/bin/python3 -c "import urllib.request as u; u.urlretrieve('https://raw.githubusercontent.com/Toki-bio/SINEplot/23234c4f3523d415e8616242f91d0d3b265aeec9/SINEplot.py', '/opt/conda/bin/SINEplot.py')" \
 && echo "fafe3e1cee97870dcb21c68fec4c15c4d0976017df210e74be7cbda8bb325958  /opt/conda/bin/SINEplot.py" | sha256sum -c - \
 && chmod 755 /opt/conda/bin/SINEplot.py

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
