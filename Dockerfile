FROM python:3.10.12

LABEL maintainer "Arianna Smith <arianna.smith@state.co.us>"
LABEL registry "Docker Hub"
LABEL namespace "ariannaesmith"
LABEL repository "cdphe_h5_influenza"
LABEL source_code "https://github.com/CDPHE-bioinformatics/CDPHE-H5-influenza"

ENV PATH="/usr/src/app:$PATH"
ENV DOCKER_VERSION='v1.2.0'
ENV APPDIR=/usr/src/app

WORKDIR $APPDIR

COPY src/scripts/* ./
COPY references/* ./references/
COPY requirements.txt ./

RUN pip3 install -r requirements.txt

CMD ["python3", "./test.py"]