# Use an official base image (e.g., Debian)
FROM debian:10

# Set the working directory
WORKDIR /app

# Install dependencies (adjust based on your application)

RUN apt-get update && \
    apt-get install -y make libblas-dev liblapacke-dev gcc vim fish

# Copy the local code to the container
COPY . .

# sudo docker build -t your_image_name:tag .
# sudo docker run -it --rm -v /home/remyy/Desktop/M1-CHPS-tp-cn-2023:/app your_image_name:tag
