#!/bin/bash
# Cloud-init user-data script for CosMx EC2 instances.
# Installs DCV, UV, Python environments, and (in analytics mode) R/RStudio.
# DCV runs a CONSOLE session on a real Xorg :0 (GDM3-managed, dummy software
# driver) so multiple Napari OpenGL windows render — unlike virtual Xdcv
# sessions, which could not. Software rendering; GPU acceleration is a
# follow-up (Phase 2).
# EC2_MODE controls what gets installed: "analytics" (default) or "napari".
# Run by create_ami.py or start_ec2.py --raw as user-data.
set -euxo pipefail

exec > >(tee /var/log/ami-setup.log) 2>&1
echo "=== AMI setup started at $(date -u) ==="

export DEBIAN_FRONTEND=noninteractive
EC2_MODE="${EC2_MODE:-analytics}"
echo "EC2_MODE=$EC2_MODE"

# ── System packages ──────────────────────────────────────────────────────
apt-get update && apt-get upgrade -y
apt-get install -y --no-install-recommends \
    curl wget git unzip tmux \
    build-essential \
    software-properties-common \
    libcurl4-openssl-dev libssl-dev libxml2-dev

# ── AWS CLI v2 (not available via apt on Ubuntu 24.04) ────────────────
curl -sSL "https://awscli.amazonaws.com/awscli-exe-linux-$(uname -m).zip" -o /tmp/awscliv2.zip
unzip -q /tmp/awscliv2.zip -d /tmp
/tmp/aws/install
rm -rf /tmp/aws /tmp/awscliv2.zip

# ── Desktop environment (for a DCV console session) ─────────────────────
# We run a real Xorg on :0 (dummy software driver) managed by GDM3, and DCV
# attaches a CONSOLE session to it. This mirrors the pre-e4926e7 setup that
# handled multiple Napari windows, but uses GDM3 (AWS's recommended display
# manager for Ubuntu 24.04) instead of LightDM, which is what broke before.
# gdm3 + xserver-xorg-video-dummy give the :0 Xorg; x11-xserver-utils = xrandr.
apt-get install -y --no-install-recommends \
    xfce4 xfce4-terminal \
    gdm3 \
    xserver-xorg-core xserver-xorg-video-dummy x11-xserver-utils \
    xfonts-base \
    desktop-file-utils \
    mesa-utils \
    dbus-x11 \
    firefox \
    libgles2 libegl1 libgl1-mesa-dri \
    libxcb-cursor0 libxcb-xinerama0 libxcb-randr0 libxcb-shape0 \
    libxcb-xfixes0 libxcb-icccm4 libxcb-image0 libxcb-keysyms1 \
    libxcb-render-util0 libxkbcommon-x11-0 \
    at-spi2-core xdg-desktop-portal xdg-desktop-portal-gtk

if [ "$EC2_MODE" != "napari" ]; then
# ── R from Ubuntu repos (pre-built binaries, fast) ──────────────────────
apt-get install -y --no-install-recommends \
    r-base \
    gfortran \
    libgmp-dev libmpfr-dev libarmadillo-dev libglpk-dev libnlopt-dev libgsl-dev \
    r-cran-ggplot2 \
    r-cran-dplyr \
    r-cran-reticulate \
    r-cran-remotes

# InSituType CRAN dependencies
Rscript -e 'install.packages(c("fastglm", "irlba", "mclust", "umap", "uwot", "RcppArmadillo", "data.table", "Matrix"), repos = "https://cloud.r-project.org")'

# InSituType Bioconductor dependencies
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org"); BiocManager::install(c("SingleCellExperiment", "SummarizedExperiment", "sparseMatrixStats", "spatstat.geom"), ask = FALSE)'

# InSituType itself (GitHub only)
Rscript -e 'remotes::install_github("Nanostring-Biostats/insitutype")'

# ── RStudio Server ──────────────────────────────────────────────────────
RSTUDIO_VERSION="2024.12.1-563"
wget -q "https://download2.rstudio.org/server/jammy/amd64/rstudio-server-${RSTUDIO_VERSION}-amd64.deb"
apt-get install -y "./rstudio-server-${RSTUDIO_VERSION}-amd64.deb" || {
    apt-get --fix-broken install -y
    dpkg -i "rstudio-server-${RSTUDIO_VERSION}-amd64.deb"
}
rm -f rstudio-server-*.deb
systemctl enable rstudio-server
fi

# ── NICE DCV (console session on a real Xorg :0, software rendering) ─────
cd /tmp
OS_VERSION=$(. /etc/os-release; echo "$VERSION_ID" | sed 's/\\.//g')
ARCH=$(arch)
DCV_TGZ="nice-dcv-ubuntu${OS_VERSION}-${ARCH}.tgz"

wget -q "https://d1uj6qtbmh3dt5.cloudfront.net/${DCV_TGZ}" || {
    echo "DCV package for Ubuntu ${OS_VERSION} not found, trying 2204 fallback..."
    DCV_TGZ="nice-dcv-ubuntu2204-${ARCH}.tgz"
    wget -q "https://d1uj6qtbmh3dt5.cloudfront.net/${DCV_TGZ}"
}

tar xzf "$DCV_TGZ"
cd nice-dcv-*-"${ARCH}"
# Console sessions use the system Xorg, so nice-xdcv (DCV's own X server for
# virtual sessions) is not needed.
apt-get install -y ./nice-dcv-server_*.deb ./nice-dcv-web-viewer_*.deb
cd /tmp && rm -rf nice-dcv-*-"${ARCH}" "$DCV_TGZ"

# Configure DCV to auto-create a CONSOLE session owned by ubuntu, attached to
# the Xorg running on :0 (started by GDM3 below). A console session uses a
# real Xorg + window manager, which — unlike a virtual Xdcv session — happily
# drives several Napari OpenGL windows at once (the workflow that regressed
# when we moved to virtual sessions).
cat > /etc/dcv/dcv.conf << 'DCVCONF'
[display]
target-fps = 30

[session-management]
create-session = true

[session-management/automatic-console-session]
owner = "ubuntu"
storage-root = "%home%"

[connectivity]
web-port = 8443
DCVCONF

systemctl enable dcvserver

# ── Xorg on :0 via GDM3 (no LightDM — that's what broke on Ubuntu 20.x+) ──
# Software framebuffer via the dummy driver at a large virtual resolution.
cat > /etc/X11/xorg.conf << 'XORGCONF'
Section "Device"
    Identifier "DummyDevice"
    Driver "dummy"
    Option "UseEDID" "false"
    VideoRam 512000
EndSection

Section "Monitor"
    Identifier "DummyMonitor"
    HorizSync   5.0 - 1000.0
    VertRefresh 5.0 - 200.0
    Option "ReducedBlanking"
EndSection

Section "Screen"
    Identifier "DummyScreen"
    Device "DummyDevice"
    Monitor "DummyMonitor"
    DefaultDepth 24
    SubSection "Display"
        Viewport 0 0
        Depth 24
        Virtual 3840 2160
    EndSubSection
EndSection
XORGCONF

# GDM3: force Xorg (DCV can't use Wayland), and auto-login ubuntu into an Xfce
# session so a desktop is running on :0 for DCV's console session to attach to.
echo "/usr/sbin/gdm3" > /etc/X11/default-display-manager
cat > /etc/gdm3/custom.conf << 'GDMCONF'
[daemon]
WaylandEnable=false
AutomaticLoginEnable=true
AutomaticLogin=ubuntu
GDMCONF
# Make the auto-login session Xfce.
mkdir -p /var/lib/AccountsService/users
cat > /var/lib/AccountsService/users/ubuntu << 'ACCT'
[User]
Session=xfce
XSession=xfce
SystemAccount=false
ACCT
systemctl enable gdm3
systemctl set-default graphical.target

# ── rclone (fast parallel S3 transfers) ─────────────────────────────────
curl -sSL https://rclone.org/install.sh | bash

# ── UV (Python package manager) ─────────────────────────────────────────
curl -LsSf https://astral.sh/uv/install.sh | sh
cp /root/.local/bin/uv /usr/local/bin/uv
cp /root/.local/bin/uvx /usr/local/bin/uvx

# ── Clone repo and install Python environments ──────────────────────────
REPO_DIR="/opt/cosmx-utilities"
GIT_BRANCH="${GIT_BRANCH:-main}"
git clone -b "$GIT_BRANCH" https://github.com/UW-BRaIN-lab/cosmx-utilities.git "$REPO_DIR"
chown -R ubuntu:ubuntu "$REPO_DIR"

# Main workspace (napari-cosmx-fork + pipeline tools) — requires Python <3.11
cd "$REPO_DIR"
sudo -u ubuntu uv python install 3.10
sudo -u ubuntu uv sync --python 3.10 --extra gui

if [ "$EC2_MODE" != "napari" ]; then
# Analytics environment (Jupyter + Polars) — uses latest Python
cd "$REPO_DIR/ec2/analytics"
sudo -u ubuntu uv sync

# ── Default mount point for data ────────────────────────────────────────
mkdir -p /mnt/cosmx
chown ubuntu:ubuntu /mnt/cosmx
fi

# ── Bring up the graphical target: GDM3 starts Xorg :0, autologin ubuntu, ─
#    then DCV auto-creates the console session against it. ─────────────────
systemctl start dcvserver
# Switch to graphical.target so GDM3 starts the :0 Xorg + Xfce session now
# (idempotent on future boots via the enabled units above).
systemctl isolate graphical.target || echo "WARNING: could not isolate graphical.target"

# ── Sentinel: signal setup completion ────────────────────────────────────
touch /var/lib/cloud/instance/ami-setup-complete
echo "=== AMI_SETUP_COMPLETE ==="
echo "=== AMI setup finished at $(date -u) ==="
