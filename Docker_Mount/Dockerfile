# Start from the official FEniCS base image
FROM dolfinx/lab:stable

# Upgrade pip and install ldrb
RUN python3 -m pip install --upgrade pip \
    && pip install fenicsx-ldrb cardiac-geometriesx numpy pyvista notebook jupyterlab meshio trimesh pymeshfix

# Install Xvfb and other utilities
RUN apt-get update && apt-get install -y xvfb x11-utils \
    && rm -rf /var/lib/apt/lists/*


# Set the working directory in the container
WORKDIR /home/fenics/shared

# Copy tester script into the container
# COPY tester.py .

# Default command to run the script
# CMD ["python3", "tester.py"]

ENV XDG_RUNTIME_DIR=/tmp/runtime-fenics
RUN mkdir -p /tmp/runtime-fenics

