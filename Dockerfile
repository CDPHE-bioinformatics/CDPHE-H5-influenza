FROM ubuntu:22.04

LABEL maintainer "Arianna Smith <arianna.smith@state.co.us>"

WORKDIR /app

COPY . .

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
  && apt-get install -y --no-install-recommends git-all

CMD ["/bin/bash"]
