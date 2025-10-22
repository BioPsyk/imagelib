# Beagle 5.5 Docker Image
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install Java runtime and dependencies
RUN apt-get update && apt-get install -y \
    openjdk-21-jre-headless \
    wget \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Download Beagle 5.5 JAR file
WORKDIR /opt/beagle
RUN wget -q https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar && \
    chmod 644 beagle.27Feb25.75f.jar

# Create wrapper script for easy execution
RUN echo '#!/bin/bash\n\
if [ $# -eq 0 ]; then\n\
    java -jar /opt/beagle/beagle.27Feb25.75f.jar\n\
else\n\
    java -Xmx4g -jar /opt/beagle/beagle.27Feb25.75f.jar "$@"\n\
fi' > /usr/local/bin/beagle && \
    chmod +x /usr/local/bin/beagle

# Copy VERSION file and create version information
COPY VERSION /etc/beagle-version-tag
RUN java -jar /opt/beagle/beagle.27Feb25.75f.jar 2>&1 | head -5 > /etc/beagle-version || \
    echo "Beagle 5.5 (27Feb25.75f)" > /etc/beagle-version

# Set working directory for data
WORKDIR /data

# Create a non-root user for security
RUN groupadd -r beagle && useradd -r -g beagle beagle && \
    chown -R beagle:beagle /data
USER beagle

# Default command shows version and help
CMD ["beagle"]
