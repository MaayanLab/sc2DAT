FROM python:3.8

ENV DEBIAN_FRONTEND "noninteractive"
ENV TZ "America/New_York"

RUN set -x \
  && echo "Preparing system..." \
  && apt-get -y update \
  && apt-get -y install \
    curl \
    fuse \
    git \
    nginx \
  && rm -rf /var/lib/apt/lists/* \
  && pip3 install --no-cache-dir --upgrade pip

ENV PLAYWRIGHT_BROWSERS_PATH /ms-playwright

RUN set -x \
  && echo "Installing jupyter kernel..." \
  && pip3 install --no-cache-dir ipython_genutils ipykernel playwright \
  && python3 -m ipykernel install \
  && python3 -m playwright install chromium

ADD setup.R /app/setup.R
RUN set -x \
  && echo "Installing R..." \
  && apt-get -y clean \
  && apt-get -y update \
  && apt-get -y install r-base \
  && rm -rf /var/lib/apt/lists/*

ADD deps.txt /app/deps.txt
RUN set -x \
  && echo "Installing system dependencies from deps.txt..." \
  && apt-get -y update \
  && apt-get -y install $(grep -v '^#' /app/deps.txt) \
  && rm -rf /var/lib/apt/lists/* \
  && rm /app/deps.txt

ADD setup.R /app/setup.R
RUN set -x \
  && echo "Setting up R with setup.R..." \
  && R -e "source('/app/setup.R')" \
  && rm /app/setup.R

ADD requirements.txt /app/requirements.txt
RUN set -x \
  && echo "Installing python dependencies from requirements.txt..." \
  && pip3 install --no-cache-dir -r /app/requirements.txt \
  && rm /app/requirements.txt

ARG appyter_version=appyter[production]@git+https://github.com/Maayanlab/appyter
RUN set -x \
  && echo "Installing appyter..." \
  && pip3 install --no-cache-dir --upgrade ${appyter_version}

RUN set -x \
  && echo "Preparing user..." \
  && useradd -ms /bin/bash -d /app app \
  && groupadd fuse \
  && adduser app fuse \
  && mkdir -p /app /app/data /data \
  && chown -R app:app /app /data \
  && chmod og+rwx -R /var/lib/nginx /var/log/nginx

USER app
WORKDIR /app
EXPOSE 5000
VOLUME /app/data

ENV APPYTER_PREFIX "/"
ENV APPYTER_HOST "0.0.0.0"
ENV APPYTER_PORT "5000"
ENV APPYTER_DEBUG "false"
ENV APPYTER_IPYNB "main.ipynb"
ENV APPYTER_PROFILE "bootstrap"
ENV APPYTER_EXTRAS '["toc", "toggle-code", "hide-code"]'

ENV PATH "/app:$PATH"
ENV PYTHONPATH "/app:$PYTHONPATH"

COPY --chown=app:app . /app

CMD [ "appyter", "flask-app" ]
