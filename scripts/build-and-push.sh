#!/bin/bash

set -e

# Configuration
REGISTRY="biopsyk"  # Docker Hub organization
PLATFORMS="linux/amd64,linux/arm64"

# Function to build and push a single image
build_and_push_image() {
    local image_name=$1
    local image_path="images/${image_name}"
    
    if [ ! -d "$image_path" ]; then
        echo "Error: Image directory $image_path not found"
        exit 1
    fi
    
    echo "Building and pushing ${image_name}..."
    
    # Set up BuildX builder if it doesn't exist
    if ! docker buildx inspect multiarch-builder >/dev/null 2>&1; then
        docker buildx create --name multiarch-builder --driver docker-container --bootstrap
    fi
    docker buildx use multiarch-builder
    
    # Build and push the multi-architecture image
    docker buildx build \
        --platform ${PLATFORMS} \
        --tag ${REGISTRY}/${image_name}:latest \
        --push \
        ${image_path}
        
    echo "Successfully built and pushed ${image_name}"
}

# Main script
if [ $# -eq 0 ]; then
    echo "Usage: $0 <image_name|all>"
    echo "Available images:"
    ls images/
    exit 1
fi

# Check if docker is logged in
if ! docker info >/dev/null 2>&1; then
    echo "Error: Docker daemon is not running"
    exit 1
fi

# Try to verify login by attempting to pull a public image
if ! docker pull hello-world:latest >/dev/null 2>&1; then
    echo "Error: Docker daemon is not responding properly"
    exit 1
fi

# Check if we can access the registry
if ! docker pull ${REGISTRY}/plink:latest >/dev/null 2>&1; then
    if ! docker pull hello-world:latest >/dev/null 2>&1; then
        echo "Error: Cannot access Docker Hub. Please run 'docker login' first."
        exit 1
    fi
fi

# Build specific image or all images
if [ "$1" = "all" ]; then
    for image in images/*; do
        if [ -d "$image" ]; then
            build_and_push_image $(basename "$image")
        fi
    done
else
    build_and_push_image "$1"
fi 