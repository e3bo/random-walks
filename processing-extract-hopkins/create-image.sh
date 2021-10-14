#!/bin/bash

set -e

# The name of our algorithm
algorithm_name=cssetransform

#set -e # stop if anything fails
account=$(aws sts get-caller-identity --query Account --output text)

# Get the region defined in the current configuration (default to us-east-1 if none defined)
region=$(aws configure get region)
region=${region:-us-east-1}

if [ "$region" = "cn-north-1" ] || [ "$region" = "cn-northwest-1" ]; then domain="amazonaws.com.cn"; 
else domain="amazonaws.com"; fi

registryname="${account}.dkr.ecr.${region}.${domain}"
fullname="${registryname}/${algorithm_name}:latest"

# If the repository doesn't exist in ECR, create it.
aws ecr describe-repositories --repository-names "${algorithm_name}" > /dev/null 2>&1

if [ $? -ne 0 ]
then
    aws ecr create-repository --repository-name "${algorithm_name}" > /dev/null
fi

aws ecr get-login-password --region ${region} | docker login --username AWS --password-stdin ${registryname}

# Build the docker image locally with the image name and then push it to ECR
# with the full name.
docker build  -t ${algorithm_name} .
docker tag ${algorithm_name} ${fullname}

docker push ${fullname}
