#!/usr/bin/env bash
set -euo pipefail

# Helper script to tag and push a release for the current version
# Checks that no tag or release already exists before proceeding
#
# Usage: bash scripts/release.sh

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

info() { echo -e "${GREEN}[INFO]${NC} $1"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

# Check we're on main branch
CURRENT_BRANCH=$(git branch --show-current)
if [[ "$CURRENT_BRANCH" != "main" ]]; then
    warn "Not on main branch (currently on '${CURRENT_BRANCH}')"
    read -p "Switch to main and pull latest? [y/N] " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        info "Switching to main..."
        git checkout main
        info "Pulling latest changes..."
        git pull origin main
    else
        error "Releases should be created from the main branch"
    fi
fi

# Ensure working directory is clean
if [[ -n $(git status --porcelain) ]]; then
    warn "Working directory has uncommitted changes"
    git status --short
    read -p "Continue anyway? [y/N] " -n 1 -r
    echo ""
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        error "Aborted - please commit or stash changes first"
    fi
fi

# Get current version from bumpver
VERSION=$(uv run bumpver show 2>/dev/null | grep "Current Version:" | awk '{print $3}')
if [[ -z "$VERSION" ]]; then
    error "Could not determine current version from bumpver"
fi

TAG="v${VERSION}"
info "Current version: ${VERSION}"
info "Tag to create: ${TAG}"

# Check if local tag exists
if git tag -l "$TAG" | grep -q "^${TAG}$"; then
    error "Local tag '${TAG}' already exists"
fi
info "No local tag '${TAG}' found"

# Check if remote tag exists
if git ls-remote --tags origin "$TAG" 2>/dev/null | grep -q "$TAG"; then
    error "Remote tag '${TAG}' already exists on origin"
fi
info "No remote tag '${TAG}' found"

# Check if GitHub release exists (requires gh CLI)
if command -v gh &> /dev/null; then
    if gh release view "$TAG" &> /dev/null; then
        error "GitHub release '${TAG}' already exists"
    fi
    info "No GitHub release '${TAG}' found"
else
    warn "gh CLI not found, skipping GitHub release check"
fi

# Confirm with user
echo ""
echo -e "${YELLOW}Ready to create and push tag '${TAG}'${NC}"
read -p "Proceed? [y/N] " -n 1 -r
echo ""

if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    info "Aborted by user"
    exit 0
fi

# Create and push the tag
info "Creating tag '${TAG}'..."
git tag -a "$TAG" -m "Release ${VERSION}"

info "Pushing tag '${TAG}' to origin..."
git push origin "$TAG"

echo ""
info "Successfully created and pushed tag '${TAG}'"
info "GitHub Actions will now build and create the release"
