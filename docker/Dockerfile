FROM docker.io/rocker/geospatial:4.0.0-ubuntu18.04

# BEGIN copying from cml Dockerfile

RUN apt-get update && \
    apt-get install -y --no-install-recommends --fix-missing \
        build-essential \
        apt-utils \
        apt-transport-https \
        ca-certificates \
        software-properties-common \
        pkg-config \
        curl \
        wget \
        unzip \
        gpg-agent \
        locales && \
    add-apt-repository universe -y && \
    add-apt-repository ppa:git-core/ppa -y && \
    apt-get install -y git && \
    curl -sL https://deb.nodesource.com/setup_12.x | bash && \
    apt-get install -y nodejs  && \
    locale-gen en_US.UTF-8 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV LANG='en_US.UTF-8' LANGUAGE='en_US:en' LC_ALL='en_US.UTF-8'

# Install update pip and nodejs, install dvc and cml
RUN wget https://dvc.org/deb/dvc.list -O /etc/apt/sources.list.d/dvc.list && \
    apt update && \
    apt -y install dvc && \
    apt -y install build-essential libcairo2-dev libpango1.0-dev libjpeg-dev libgif-dev librsvg2-dev && \
    npm config set user 0 && \
    npm config set unsafe-perm true && \
    npm install -g @dvcorg/cml && \
    npm install -g vega-cli && \
    npm install -g vega-lite && \
    apt-get install -y libfontconfig-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# DOCKER, OUR CUSTOM DOCKER MACHINE FORK, GITLAB RUNNER AND GITHUB RUNNER
# OUR DOCKER_MACHINE FORK SUPPORTS GPU ACCELARATOR
ENV RUNNER_PATH=/home/runner
ENV RUNNER_ALLOW_RUNASROOT=1

# FIXES Python setup Action
RUN mkdir /__t
RUN ln -s /__t /opt/hostedtoolcache
ENV AGENT_TOOLSDIRECTORY=/opt/hostedtoolcache

RUN mkdir ${RUNNER_PATH}
WORKDIR ${RUNNER_PATH}

RUN wget -O /usr/local/bin/gitlab-runner https://gitlab-runner-downloads.s3.amazonaws.com/latest/binaries/gitlab-runner-linux-amd64 && \
    chmod +x /usr/local/bin/gitlab-runner && \
    gitlab-runner install --user=root --working-directory=${RUNNER_PATH} && \
    wget https://github.com/actions/runner/releases/download/v2.267.1/actions-runner-linux-x64-2.267.1.tar.gz && \
    tar xzf actions-runner-linux-x64-2.267.1.tar.gz && \
    ./bin/installdependencies.sh && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

CMD ["cml-cloud-runner-entrypoint"]

# END copying from cml Dockerfile

RUN apt-get update && \
    apt-get install -y --no-install-recommends --fix-missing \
        subversion

RUN install2.r --skipinstalled --error \
ggforce \
JuliaCall \
KScorrect \
lbfgs \
lme4 \
padr \
tictoc \
verification

WORKDIR /opt
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.3-linux-x86_64.tar.gz
RUN tar zxvf julia-1.5.3-linux-x86_64.tar.gz
ENV PATH="${PATH}:/opt/julia-1.5.3/bin"
USER rstudio
RUN /opt/julia-1.5.3/bin/julia -e 'using Pkg; Pkg.add(name = "MbedTLS", version = "0.7.0"); Pkg.add.(["CSVFiles", "DataFrames", "Distributions", "ForwardDiff", "LinearAlgebra", "Optim"])'
USER root
WORKDIR ${RUNNER_PATH}
RUN apt-get update && apt-get install --only-upgrade dvc
RUN installGithub.r reichlab/zoltr epiforecasts/scoringutils reichlab/covidHubUtils
