FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Install only what we need and keep image small.
# Pin numpy to 1.x to avoid pyBigWig import errors under NumPy 2.x.
RUN pip install --no-cache-dir pysam==0.22.1

RUN mkdir -p /opt2 /data2
COPY Dockerfile /opt2/Dockerfile
RUN chmod a+rX /opt2/Dockerfile
WORKDIR /data2

ENTRYPOINT ["/bin/bash"]
