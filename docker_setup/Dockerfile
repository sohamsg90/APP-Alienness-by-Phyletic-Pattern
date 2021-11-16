FROM ubuntu:16.04
WORKDIR /home
RUN apt-get update && \
    apt-get install build-essential -y && \
    apt-get install libimage-magick-perl -y && \
    apt-get install wget -y && \
    apt-get install curl -y && \
    apt-get install cpanminus -y && \
    apt-get install bash-builtins -y && \
    apt-get install dirmngr gnupg apt-transport-https ca-certificates software-properties-common -y && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' && \
    apt-get install r-base -y && \
    apt-get install r-base-dev -y && \
    apt-get install git -y && \
    apt-get install hmmer -y && \
    apt-get install default-jre -y

ENV HOME_PATH=/home

COPY files_to_export/* /home/
COPY APP_scripts/ /home/APP_scripts/
COPY reserve/ /home/reserve/
RUN sh run.sh


