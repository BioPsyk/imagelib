FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# Install a full Perl interpreter (not just perl-base)
RUN apt-get update && apt-get install -y \
    perl \
    && rm -rf /var/lib/apt/lists/*

# Copy VERSION file and record installed Perl version
COPY VERSION /etc/perl-version-tag
RUN echo "Perl $(perl -e 'print $^V') - Image $(cat /etc/perl-version-tag)" > /etc/perl-version

WORKDIR /data

# Create a non-root user for security
RUN groupadd -r perl && useradd -r -g perl perl
USER perl

CMD ["perl", "-v"]
