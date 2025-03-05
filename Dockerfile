FROM python:3.10.12

LABEL maintainer "Arianna Smith <arianna.smith@state.co.us>"

ENV PATH="/usr/src/app:$PATH"
ENV DOCKER_VERSION='v0.1.0'
ENV APPDIR=/usr/src/app

WORKDIR $APPDIR

COPY src/scripts/* ./
COPY references/* ./references/
COPY requirements.txt ./

RUN pip3 install -r requirements.txt

CMD ["python3", "./test.py"]