FROM rocker/tidyverse:4.1.2

RUN apt-get update -qq && \
    # Clean up after apt.
    rm -rf /var/lib/apt/lists/*

RUN install2.r -e \
      ggpmisc \
      patchwork
