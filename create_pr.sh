#!/bin/bash

# Script to create a pull request for the Documenter.jl documentation

echo "This script will help you create a pull request for the documentation."
echo ""
echo "Prerequisites:"
echo "1. You need to have a GitHub repository for CurveFit.jl"
echo "2. You need to have the GitHub CLI (gh) installed"
echo "3. You need to be authenticated with gh (run 'gh auth login' if not)"
echo ""

read -p "Have you set up your GitHub repository? (y/n): " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Please create a GitHub repository first, then run this script again."
    exit 1
fi

echo ""
read -p "Enter your GitHub repository URL (e.g., https://github.com/username/CurveFit.jl): " REPO_URL

# Extract owner and repo name from URL
REPO_NAME=$(basename "$REPO_URL" .git)
OWNER=$(basename $(dirname "$REPO_URL"))

# Add remote if not already added
if ! git remote | grep -q origin; then
    echo "Adding remote origin..."
    git remote add origin "$REPO_URL"
fi

# Push the main branch first if needed
echo "Pushing main branch..."
git push -u origin master:main 2>/dev/null || git push -u origin main

# Push the documentation branch
echo "Pushing documentation branch..."
git push -u origin add-documenter-docs

# Create the pull request using gh CLI
echo ""
echo "Creating pull request..."
gh pr create \
    --title "Add Documenter.jl documentation for CurveFit.jl" \
    --body-file PR_DESCRIPTION.md \
    --base main \
    --head add-documenter-docs \
    --repo "$OWNER/$REPO_NAME"

echo ""
echo "Pull request created successfully!"
echo "You can view it at: https://github.com/$OWNER/$REPO_NAME/pulls"