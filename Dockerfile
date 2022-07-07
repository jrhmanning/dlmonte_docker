FROM python:3.8-buster

RUN apt-get update && \
    cat /etc/apt/sources.list && \
    apt-get -y install make && \
    make -v && \
    apt-get -y install wget && \
	apt-get -y install gfortran && \
	gfortran -v && \
	python3 --version && \
	pip --version && \
	pip install scipy && \
	pip install numpy && \
	pip install pandas && \
	pip install PyYAML && \
	pip install dlmontepython && \
	pip install ase && \
	pip install matplotlib && \
	pip freeze
	
RUN wget -q https://gitlab.com/dl_monte/DL_MONTE-2/-/archive/master/DL_MONTE-2-master.zip && \
    unzip DL_MONTE-2-master.zip && \
	cd DL_MONTE-2-master && \
	bash build: SRL dir gfortranrepro -f bin/DLMONTE-SRL.X && \
	mv bin/DLMONTE-SRL.X /run/DLMONTE-SRL.X

COPY interface/*.cif /run
COPY scripts/*.py /run

WORKDIR /run

CMD ["python", "isotherm_runner.py"]
