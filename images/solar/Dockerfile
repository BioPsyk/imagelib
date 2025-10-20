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

# Download and install SOLAR-Eclipse 9.0.1 static version
WORKDIR /tmp
RUN wget -q https://www.nitrc.org/frs/download.php/12509/solar-eclipse-9.0.1-static-Linux.zip && \
    unzip -q solar-eclipse-9.0.1-static-Linux.zip && \
    mkdir -p /opt/solar-eclipse && \
    mv solar-eclipse-9.0.1-static-Linux/* /opt/solar-eclipse/ && \
    rm -rf /tmp/*

# Add SOLAR to PATH
ENV PATH="/opt/solar-eclipse/bin:${PATH}"

# Copy VERSION file and create version information
COPY VERSION /etc/solar-version-tag
RUN echo "SOLAR Eclipse 9.0.1 (Static Linux)" > /etc/solar-version

# Set working directory for data
WORKDIR /data

# Create a non-root user for security
RUN groupadd -r solar && useradd -r -g solar solar && \
    chown -R solar:solar /data
USER solar

# Default command shows help
CMD ["solar"]
