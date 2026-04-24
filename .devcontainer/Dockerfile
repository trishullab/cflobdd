FROM node:20

ARG TZ
ENV TZ="$TZ"

ARG CLAUDE_CODE_VERSION=latest

# Install basic development tools and iptables/ipset
RUN apt-get update && apt-get install -y --no-install-recommends \
  less \
  git \
  procps \
  sudo \
  fzf \
  zsh \
  man-db \
  unzip \
  gnupg2 \
  gh \
  iptables \
  ipset \
  iproute2 \
  dnsutils \
  aggregate \
  jq \
  nano \
  vim \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

# CFLOBDD project-specific dependencies.
# Kept in a separate RUN layer so local changes here don't invalidate
# the cache of the upstream Anthropic-reference apt block above.
RUN apt-get update && apt-get install -y --no-install-recommends \
  build-essential \
  gdb \
  cmake \
  lld \
  libboost-all-dev \
  flex \
  python3-pip \
  python3-dev \
  python3-matplotlib \
  texlive-latex-recommended \
  texlive-latex-extra \
  latexmk \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

# Intel oneAPI VTune Profiler (Linux).
# Only Hotspots (user-mode sampling) is supported in this container,
# which does not need PMU / hardware-event sampling. See devcontainer.json
# runArgs for the SYS_PTRACE cap and seccomp setting required at runtime.
RUN wget -qO- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
    | gpg --dearmor -o /usr/share/keyrings/oneapi-archive-keyring.gpg && \
  echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
    > /etc/apt/sources.list.d/oneAPI.list && \
  apt-get update && apt-get install -y --no-install-recommends \
    intel-oneapi-vtune \
  && apt-get clean && rm -rf /var/lib/apt/lists/*

# Put the VTune CLI on PATH without sourcing the full oneAPI setvars.sh
# (which is verbose and sets many unrelated vars).
ENV PATH="/opt/intel/oneapi/vtune/latest/bin64:${PATH}"

# Ensure default node user has access to /usr/local/share
RUN mkdir -p /usr/local/share/npm-global && \
  chown -R node:node /usr/local/share

ARG USERNAME=node

# Persist bash history.
RUN SNIPPET="export PROMPT_COMMAND='history -a' && export HISTFILE=/commandhistory/.bash_history" \
  && mkdir /commandhistory \
  && touch /commandhistory/.bash_history \
  && chown -R $USERNAME /commandhistory

# Set `DEVCONTAINER` environment variable to help with orientation
ENV DEVCONTAINER=true

# Create workspace and config directories and set permissions.
# The workspaceFolder in devcontainer.json mounts the project at
# /mnt/c/Users/reps/Documents/Software/CFLOBDD-github/cflobdd to match
# the host path, so create that directory hierarchy with node ownership.
RUN mkdir -p /mnt/c/Users/reps/Documents/Software/CFLOBDD-github/cflobdd /home/node/.claude && \
  chown -R node:node /mnt/c /home/node/.claude

WORKDIR /mnt/c/Users/reps/Documents/Software/CFLOBDD-github/cflobdd

ARG GIT_DELTA_VERSION=0.18.2
RUN ARCH=$(dpkg --print-architecture) && \
  wget "https://github.com/dandavison/delta/releases/download/${GIT_DELTA_VERSION}/git-delta_${GIT_DELTA_VERSION}_${ARCH}.deb" && \
  sudo dpkg -i "git-delta_${GIT_DELTA_VERSION}_${ARCH}.deb" && \
  rm "git-delta_${GIT_DELTA_VERSION}_${ARCH}.deb"

# Set up non-root user
USER node

# Install global packages
ENV NPM_CONFIG_PREFIX=/usr/local/share/npm-global
ENV PATH=$PATH:/usr/local/share/npm-global/bin

# Set the default shell to zsh rather than sh
ENV SHELL=/bin/zsh

# Set the default editor and visual
ENV EDITOR=nano
ENV VISUAL=nano

# Default powerline10k theme
ARG ZSH_IN_DOCKER_VERSION=1.2.0
RUN sh -c "$(wget -O- https://github.com/deluan/zsh-in-docker/releases/download/v${ZSH_IN_DOCKER_VERSION}/zsh-in-docker.sh)" -- \
  -p git \
  -p fzf \
  -a "source /usr/share/doc/fzf/examples/key-bindings.zsh" \
  -a "source /usr/share/doc/fzf/examples/completion.zsh" \
  -a "export PROMPT_COMMAND='history -a' && export HISTFILE=/commandhistory/.bash_history" \
  -x

# Install Claude
RUN npm install -g @anthropic-ai/claude-code@${CLAUDE_CODE_VERSION}


# Copy and set up firewall script
COPY init-firewall.sh /usr/local/bin/
USER root
RUN chmod +x /usr/local/bin/init-firewall.sh && \
  echo "node ALL=(root) NOPASSWD: /usr/local/bin/init-firewall.sh" > /etc/sudoers.d/node-firewall && \
  chmod 0440 /etc/sudoers.d/node-firewall
USER node
