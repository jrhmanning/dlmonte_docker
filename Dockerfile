FROM python:3.8-slim-buster
WORKDIR /run

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gfortran \
    make \
    unzip \
    wget \
  && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY Pipfile Pipfile
RUN python -m pip install --no-cache-dir --upgrade pip pipenv
RUN pipenv install --skip-lock --system --verbose && \
    pipenv --clear

# Install dlmontepython
RUN wget -q https://gitlab.com/dl_monte/dlmontepython/-/archive/master/dlmontepython-master.zip && \
    unzip dlmontepython-master.zip && \
    cd dlmontepython-master && \
    python setup.py install && \
    cd .. && \
    rm -fr dlmontepython-master.zip dlmontepython-master
	
	
# Install DLMONTE
RUN wget -q https://gitlab.com/dl_monte/DL_MONTE-2/-/archive/master/DL_MONTE-2-master.zip && \
    unzip DL_MONTE-2-master.zip && \
    cd DL_MONTE-2-master && \
      bash build: SRL dir gfortranrepro -f bin/DLMONTE-SRL.X && \
      mv bin/DLMONTE-SRL.X /usr/local/bin/DLMONTE-SRL.X && \
    cd .. && \
    rm -fr DL_MONTE-2-master.zip DL_MONTE-2-master

# Copy Python scripts to image
COPY scripts/*.py /run