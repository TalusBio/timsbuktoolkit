FROM debian:bullseye-slim

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps ca-certificates && \
    update-ca-certificates && \
    rm -rf /var/lib/apt /var/lib/dpkg /var/lib/cache /var/lib/log

WORKDIR /app

COPY target/x86_64-unknown-linux-gnu/release/timsquery /app/timsquery
COPY target/x86_64-unknown-linux-gnu/release/timsseek /app/timsseek

ENV PATH="/app:$PATH"
