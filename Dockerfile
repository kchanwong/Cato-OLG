FROM julia:1.12-bookworm

# awscli reads the job JSON from s3:// and writes the results to $OUTPUT_S3_URI.
RUN apt-get update \
 && apt-get install -y --no-install-recommends awscli \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/olg
ENV JULIA_PROJECT=/opt/olg JULIA_NUM_THREADS=auto OUTPUT_DIR=/data

# Resolve and precompile before the source, so that a change to the model does
# not rebuild the packages.
COPY Project.toml Manifest.toml ./
RUN julia -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

COPY source/ source/
COPY cached_objects.jld2 ./
COPY run_job.jl ./

RUN mkdir /data && julia run_job.jl --selftest

ENTRYPOINT ["julia", "run_job.jl"]
