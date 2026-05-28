#!/bin/bash
# Cloud-init user-data script for CosMx EC2 instances.
# Installs DCV, UV, Python environments, and (in analytics mode) R/RStudio.
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

# ── Desktop environment (lightweight, for DCV) ──────────────────────────
apt-get install -y --no-install-recommends \
    xfce4 xfce4-terminal \
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

# ── NICE DCV (software rendering, no GPU) ───────────────────────────────
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
apt-get install -y ./nice-dcv-server_*.deb ./nice-dcv-web-viewer_*.deb ./nice-xdcv_*.deb
cd /tmp && rm -rf nice-dcv-* "$DCV_TGZ"

# Configure DCV for virtual sessions (no display manager needed)
cat > /etc/dcv/dcv.conf << 'DCVCONF'
[display]
target-fps = 30

[connectivity]
web-port = 8443
DCVCONF

systemctl enable dcvserver

# Systemd service to auto-create a DCV virtual session on boot
cat > /etc/systemd/system/dcv-virtual-session.service << 'UNIT'
[Unit]
Description=Create DCV virtual session
After=dcvserver.service
Requires=dcvserver.service

[Service]
Type=oneshot
ExecStartPre=/bin/sleep 2
ExecStart=/usr/bin/dcv create-session --type virtual --owner ubuntu main
RemainAfterExit=yes

[Install]
WantedBy=multi-user.target
UNIT
systemctl enable dcv-virtual-session

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

# ── Start DCV now (services are enabled for future boots) ─────────────
systemctl start dcvserver
sleep 2
systemctl start dcv-virtual-session

# ── Sentinel: signal setup completion ────────────────────────────────────
touch /var/lib/cloud/instance/ami-setup-complete
echo "=== AMI_SETUP_COMPLETE ==="
echo "=== AMI setup finished at $(date -u) ==="
