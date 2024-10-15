FROM python:3.10.12

LABEL maintainer "Arianna Smith <arianna.smith@state.co.us>"

WORKDIR /usr/src/app

ENV PATH="$PATH:/usr/src/app"

ENV DOCKER_VERSION='v0.0.0-epsilon'

COPY src/scripts/* ./

COPY requirements.txt ./

RUN pip3 install -r requirements.txt

CMD ["python3", "./test.py"]