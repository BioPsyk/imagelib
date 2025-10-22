# SOLAR-Eclipse Docker Image
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    ca-certificates \
    tcsh \
    libgomp1 \
    libgfortran5 \
    && rm -rf /var/lib/apt/lists/*

# Download and install SOLAR-Eclipse 9.0.0 static version
WORKDIR /tmp
RUN wget -q https://www.nitrc.org/frs/download.php/12468/solar-eclipse-9.0.0-static-linux.zip && \
    unzip -q solar-eclipse-9.0.0-static-linux.zip && \
    cd solar900 && \
    ./install_solar /opt/solar-eclipse /usr/local/bin && \
    cd /tmp && \
    rm -rf /tmp/*

# Copy VERSION file and create version information
COPY VERSION /etc/solar-version-tag
RUN echo "SOLAR Eclipse 9.0.0 (Static Linux)" > /etc/solar-version

# Set working directory for data
WORKDIR /data

# Create a non-root user for security
RUN groupadd -r solar && useradd -r -g solar solar && \
    chown -R solar:solar /data
USER solar

# Default command shows help
CMD ["solar"]
