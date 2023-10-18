FROM python:3.8-slim-buster
WORKDIR /run

# Install system dependencies
RUN apt-get update && apt-get install -y \
    gfortran \
    make \
    unzip \
    wget \
	 vim \
	 nano \
	 gnuplot \
  && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY Pipfile Pipfile
RUN python -m pip install --no-cache-dir --upgrade pip pipenv
RUN pipenv install --skip-lock --system --verbose && \
    pipenv --clear

# Install dlmontepython
<<<<<<< Updated upstream
RUN wget -q https://gitlab.com/dl_monte/dlmontepython/-/archive/master/dlmontepython-master.zip && \
    unzip dlmontepython-master.zip && \
    cd dlmontepython-master && \
      python setup.py install && \
=======
RUN wget -q https://gitlab.com/dl_monte/dlmontepython/-/archive/better_fe_sweeps/dlmontepython-better_fe_sweeps.zip && \
    unzip dlmontepython-better_fe_sweeps.zip && \
    cd dlmontepython-better_fe_sweeps && \
    python setup.py install && \
>>>>>>> Stashed changes
    cd .. && \
    rm -fr dlmontepython-better_fe_sweeps.zip
	
	
# Install DLMONTE
RUN wget -q https://gitlab.com/dl_monte/DL_MONTE-2/-/archive/master/DL_MONTE-2-master.zip && \
    unzip DL_MONTE-2-master.zip && \
    cd DL_MONTE-2-master && \
      bash build.sh SRL dir gfortranrepro -f bin/DLMONTE-SRL.X && \
      mv bin/DLMONTE-SRL.X /usr/local/bin/DLMONTE-SRL.X && \
    cd .. && \
    rm -fr DL_MONTE-2-master.zip DL_MONTE-2-master

# Copy Python scripts to image
<<<<<<< Updated upstream
COPY scripts/*.py /run
=======
COPY scripts/*.py /run
COPY test_io /test_io
COPY dlmonte_gcmc_lj /run/dlmonte_gcmc_lj
COPY dlmonte_fed_lj /run/dlmonte_fed_lj
COPY debugging/* /debugging/
COPY dlmonte_tutorials/* /tutorials/
>>>>>>> Stashed changes
